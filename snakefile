# List of SRA samples used in the analysis
SAMPLES = ["SRR5660030", "SRR5660033", "SRR5660044", "SRR5660045"]

# Name of final output file
REPORT = "Ho_PipelineReport.txt"


# rule all: pipeline completes when final report exists
rule all:
    input:
        REPORT


# STEP 2A: Download HCMV genome + build index

# Download the HCMV reference genome (GCF_000845245.1)
# and extract the FASTA file.
rule download_reference:
    output:
        "reference/hcmv.fasta"
    shell:
        """
        # Create directory for reference genome
        mkdir -p reference

        # Download genome archive from NCBI using datasets CLI
        datasets download genome accession GCF_000845245.1 \
            --filename reference/hcmv.zip

        # Unzip downloaded archive into reference directory
        unzip -o reference/hcmv.zip -d reference/

        # Move FASTA file from extracted dataset structure
        # to expected pipeline location
        mv reference/ncbi_dataset/data/GCF_000845245.1/*.fna {output}
        """


# Build Bowtie2 index from reference genome FASTA
rule build_index:
    input:
        "reference/hcmv.fasta"
    output:
        "reference/hcmv_index.1.bt2"
    shell:
        """
        # Convert FASTA file into Bowtie2 index files
        # This generates multiple .bt2 index files
        bowtie2-build {input} reference/hcmv_index
        """


# STEP 2B: Map reads and retain only reads aligning to HCMV

# Align paired-end reads to reference genome.
# Only concordantly aligned read pairs are written out.
rule filter_reads:
    input:
        # Change sample_data directory to data to use full size fastq files
        r1="sample_data/{sample}_1.fastq",
        r2="sample_data/{sample}_2.fastq",
        idx="reference/hcmv_index.1.bt2"
    output:
        r1="results/filtered/{sample}_1.fastq",
        r2="results/filtered/{sample}_2.fastq",
        log="results/filtered/{sample}_bowtie_report.txt"
    shell:
        """
        # Create output directory for filtered reads
        mkdir -p results/filtered

        # Align reads to HCMV reference genome
        # --very-sensitive increases alignment sensitivity
        # --al-conc writes only read pairs that align concordantly
        # -S /dev/null discards SAM output (we only need filtered FASTQ)
        # stderr is redirected to a log file
        bowtie2 -x reference/hcmv_index \
          -1 {input.r1} -2 {input.r2} \
          --very-sensitive \
          --threads 4 \
          --al-conc results/filtered/{wildcards.sample}_%.fastq \
          -S /dev/null \
          2> {output.log}
        """


# STEP 3: SPAdes assembly (k=99)

# Perform de novo assembly using filtered reads.
rule spades:
    input:
        r1="results/filtered/{sample}_1.fastq",
        r2="results/filtered/{sample}_2.fastq"
    output:
        "results/assembly/{sample}/contigs.fasta"
    shell:
        """
        # Create assembly output directory for this sample
        mkdir -p results/assembly/{wildcards.sample}

        # Run SPAdes assembler
        # -k 99 sets k-mer size to 99 as required
        # --only-assembler skips read error correction
        spades.py -1 {input.r1} -2 {input.r2} \
            -k 99 --only-assembler \
            -o results/assembly/{wildcards.sample}
        """


# STEP 4: Assembly statistics (>1000 bp)

# Count number of contigs > 1000 bp
# and compute total length of those contigs.
rule assembly_stats:
    input:
        "results/assembly/{sample}/contigs.fasta"
    output:
        "results/assembly/{sample}/stats.txt"
    run:
        from pipeline_utils import assembly_stats

        c, total = assembly_stats(input[0], min_length=1000)

        with open(output[0], "w") as out:
            out.write(f"{c}\t{total}\n")


# STEP 5A: Extract longest contig

# Identify and write the single longest contig
# from each SPAdes assembly.
rule longest_contig:
    input:
        "results/assembly/{sample}/contigs.fasta"
    output:
        "results/blast/{sample}_longest.fasta"
    run:
        from pipeline_utils import write_longest_contig
        import os

        # Ensure BLAST results directory exists
        os.makedirs("results/blast", exist_ok=True)

        # Write longest contig to FASTA file
        write_longest_contig(input[0], output[0])


# STEP 5B: Build local BLAST database

# Create a nucleotide BLAST database
# from the HCMV reference genome.
rule build_betaherpes_db:
    input:
        "reference/hcmv.fasta"
    output:
        "betaherpes_db/betaherpes.nsq"
    shell:
        """
        # Create directory for BLAST database files
        mkdir -p betaherpes_db

        # Convert FASTA file into BLAST nucleotide database
        # Generates .nsq, .nin, .nhr index files
        makeblastdb \
          -in {input} \
          -dbtype nucl \
          -out betaherpes_db/betaherpes
        """


# STEP 5C: Run blastn

# Align longest contig against local BLAST database.
rule blast_contig:
    input:
        query="results/blast/{sample}_longest.fasta",
        db="betaherpes_db/betaherpes.nsq"
    output:
        "results/blast/{sample}_blast.tsv"
    shell:
        """
        # Run nucleotide BLAST
        # -max_hsps 1 keeps best alignment per query–subject pair
        # -max_target_seqs 10 returns up to 10 matching subjects
        # -outfmt 6 outputs tab-delimited format with specified columns
        blastn \
          -query {input.query} \
          -db betaherpes_db/betaherpes \
          -max_hsps 1 \
          -max_target_seqs 10 \
          -outfmt "6 sacc pident length qstart qend sstart send bitscore evalue stitle" \
          -out {output}
        """


# FINAL REPORT

# Combine read counts, assembly stats, and
# top 5 BLAST hits into Ho_PipelineReport.txt
rule final_report:
    input:
        idx="reference/hcmv_index.1.bt2",
        stats=expand("results/assembly/{sample}/stats.txt", sample=SAMPLES),
        blasts=expand("results/blast/{sample}_blast.tsv", sample=SAMPLES)
    output:
        REPORT
    run:
        from pipeline_utils import count_fastq_reads

        # Helper function to count paired-end read pairs
        def read_pairs(r1, r2):
            return min(count_fastq_reads(r1), count_fastq_reads(r2))

        with open(output[0], "w") as out:
            for sample in SAMPLES:
                # Change sample_data directory to data to use full size fastq files
                raw_r1 = f"sample_data/{sample}_1.fastq"
                raw_r2 = f"sample_data/{sample}_2.fastq"
                filt_r1 = f"results/filtered/{sample}_1.fastq"
                filt_r2 = f"results/filtered/{sample}_2.fastq"
                stats_file = f"results/assembly/{sample}/stats.txt"
                blast_file = f"results/blast/{sample}_blast.tsv"

                before_pairs = read_pairs(raw_r1, raw_r2)
                after_pairs = read_pairs(filt_r1, filt_r2)

                out.write(
                    f"Sample {sample} had {before_pairs} read pairs before and {after_pairs} read pairs after Bowtie2 filtering.\n"
                )

                with open(stats_file) as f:
                    c, total = f.read().strip().split("\t")

                out.write(
                    f"In the assembly of sample {sample}, there are {c} contigs > 1000 bp and {total} total bp.\n\n"
                )

                # Write BLAST header and top 5 hits
                out.write(f"{sample}:\n")
                out.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")

                with open(blast_file) as bf:
                    for i, line in enumerate(bf):
                        if i >= 5:
                            break
                        out.write(line)

                out.write("\n")