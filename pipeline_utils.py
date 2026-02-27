# Import necessary libraries
import os
from Bio import SeqIO


# STEP 2: Count FASTQ reads before and after Bowtie2 filtering
# FASTQ format uses 4 lines per read (header, sequence, separator, quality).
# Dividing total line count by 4 gives the number of reads.
def count_fastq_reads(filename):
    with open(filename) as f:
        return sum(1 for _ in f) // 4


# STEP 4: Compute assembly statistics for contigs > 1000 bp
def assembly_stats(contigs_file, min_length=1000):
    count = 0
    total_bp = 0

    # Parse contigs FASTA file using Biopython
    for record in SeqIO.parse(contigs_file, "fasta"):
        L = len(record.seq)

        # Only include contigs strictly greater than threshold
        if L > min_length:
            count += 1
            total_bp += L
    #   count = number of contigs longer than min_length
    #   total_bp = total base pairs across those contigs
    return count, total_bp


# STEP 5A: Extract the longest contig from the assembly
# The longest contig is used as the BLAST query sequence.
def write_longest_contig(contigs_fasta, out_fasta):
    longest = None
    max_len = 0

    # Iterate through all contigs in FASTA file
    for record in SeqIO.parse(contigs_fasta, "fasta"):
        length = len(record.seq)

        # Update if this contig is longer than previous maximum
        if length > max_len:
            max_len = length
            longest = record

    # Raise error if no contigs were found
    if longest is None:
        raise ValueError(f"No contigs found in {contigs_fasta}")

    # Write longest contig to output FASTA file
    SeqIO.write(longest, out_fasta, "fasta")