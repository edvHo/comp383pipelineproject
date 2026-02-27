# COMP 383 Python Pipeline Project

## Description

This project utilizes Snakemake to analyze Human Cytomegalovirus (HCMV; Human herpesvirus 5) transcriptomes collected at 2 and 6 days post-infection (dpi) from two patient donors. The HCMV reference genome (GCF_000845245.1) is downloaded and indexed with Bowtie2, and paired-end reads are aligned to filter for viral sequences. Filtered reads are assembled de novo using SPAdes (k-mer size = 99), and assembly statistics are calculated. The longest contig from each assembly is then queried using BLAST+ against a local database of Betaherpesvirinae genomes to identify the most closely related HCMV strains.

### Software
* [Snakemake](https://snakemake.readthedocs.io/en/stable/tutorial/basics.html#step-1-mapping-reads)
* [SRA toolkit](https://github.com/ncbi/sra-tools)
* [Fastq-dump](https://rnnh.github.io/bioinfo-notebook/docs/fastq-dump.html)
* [Bowtie2](https://github.com/BenLangmead/bowtie2)
* [SPAdes](https://github.com/ablab/spades)
* [BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279690/)

### Packages
* BioPython
* os

### Set Up Environment

* Clone this repository into your workspace
```
git clone https://github.com/edvHo/comp383pipelineproject.git
```
* cd into the pipeline_project directory
```
cd comp383pipelineproject
```
* Create the following directories within the cloned repository to set up input files (the sample_data directory should already be present after cloning the repo).
```
mkdir -p data
```
* This creates the directory that will contain the full size fastq files
```
mkdir -p sample_data
```
* If not already present from cloning the repository, the sample_data directory will contain trimmed fastq files

### Downloading Input Files
* Download paired-end read samples from NCBI SRA using fasterq-dump to the data directory
```
fasterq-dump SRR5660030 --split-files -O data
```
```
fasterq-dump SRR5660033 --split-files -O data
```
```
fasterq-dump SRR5660044 --split-files -O data
```
```
fasterq-dump SRR5660045 --split-files -O data
```
* Create trimmed data to use for testing to the sample_data directory
```
head -n 40000 data/SRR5660030_1.fastq > sample_data/SRR5660030_1.fastq
head -n 40000 data/SRR5660030_2.fastq > sample_data/SRR5660030_2.fastq
```
```
head -n 40000 data/SRR5660033_1.fastq > sample_data/SRR5660030_1.fastq
head -n 40000 data/SRR5660033_2.fastq > sample_data/SRR5660030_2.fastq
```
```
head -n 40000 data/SRR5660044_1.fastq > sample_data/SRR5660030_1.fastq
head -n 40000 data/SRR5660044_2.fastq > sample_data/SRR5660030_2.fastq
```
```
head -n 40000 data/SRR5660045_1.fastq > sample_data/SRR5660030_1.fastq
head -n 40000 data/SRR5660045_2.fastq > sample_data/SRR5660030_2.fastq
```

### Running the Pipeline
* To run the pipeline
```
snakemake --cores 4
```
* This will run snakefile and should output the file Ho_PipelineReport.txt, which may be checked with
```
cat Ho_PipelineReport.txt
```
