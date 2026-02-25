SAMPLES = ["SRR5660030", "SRR5660033", "SRR5660044", "SRR5660045"]
REPORT = "Ho_PipelineReport.txt"

rule all:
    input:
        REPORT


# STEP 2A: DOWNLOAD HCMV GENOME + BUILD INDEX
rule download_reference:
    output:
        "reference/hcmv.fasta"
    shell:
        """
        mkdir -p reference
        datasets download genome accession GCF_000845245.1 --filename reference/hcmv.zip
        unzip -o reference/hcmv.zip -d reference/
        mv reference/ncbi_dataset/data/GCF_000845245.1/*.fna {output}
        """

rule build_index:
    input:
        "reference/hcmv.fasta"
    output:
        "reference/hcmv_index.1.bt2"
    shell:
        """
        bowtie2-build {input} reference/hcmv_index
        """


# STEP 2B: MAP + FILTER READS
rule filter_reads:
    input:
        r1="sample_data/{sample}_1.fastq",
        r2="sample_data/{sample}_2.fastq",
        idx="reference/hcmv_index.1.bt2"
    output:
        r1="results/filtered/{sample}_1.fastq",
        r2="results/filtered/{sample}_2.fastq",
        log="results/filtered/{sample}_bowtie_report.txt"
    shell:
        """
        mkdir -p results/filtered

        bowtie2 -x reference/hcmv_index \
          -1 {input.r1} -2 {input.r2} \
          --very-sensitive \
          --threads 4 \
          --al-conc results/filtered/{wildcards.sample}_%.fastq \
          -S /dev/null \
          2> {output.log}
        """


# STEP 3 — SPAdes assembly (k=99)
rule spades:
    input:
        r1="results/filtered/{sample}_1.fastq",
        r2="results/filtered/{sample}_2.fastq"
    output:
        "results/assembly/{sample}/contigs.fasta"
    shell:
        """
        mkdir -p results/assembly/{wildcards.sample}
        spades.py -1 {input.r1} -2 {input.r2} \
            -k 99 --only-assembler \
            -o results/assembly/{wildcards.sample}
        """


# STEP 4 — Assembly stats (>1000 bp)
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


# FINAL REPORT (Steps 2–4)
rule final_report:
    input:
        idx="reference/hcmv_index.1.bt2",
        filtered=expand("results/filtered/{sample}_1.fastq", sample=SAMPLES)
                 + expand("results/filtered/{sample}_2.fastq", sample=SAMPLES),
        stats=expand("results/assembly/{sample}/stats.txt", sample=SAMPLES)
    output:
        REPORT
    run:
        from pipeline_utils import count_fastq_reads
        import os

        def read_pairs(r1, r2):
            # Paired-end: report read PAIRS; assume files match
            return min(count_fastq_reads(r1), count_fastq_reads(r2))

        with open(output[0], "w") as out:
            for sample in SAMPLES:
                # change sample_data to data for full reads
                raw_r1 = f"sample_data/{sample}_1.fastq"
                raw_r2 = f"sample_data/{sample}_2.fastq"
                filt_r1 = f"results/filtered/{sample}_1.fastq"
                filt_r2 = f"results/filtered/{sample}_2.fastq"
                stats_file = f"results/assembly/{sample}/stats.txt"

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