import os
from Bio import SeqIO

def count_fastq_reads(filename):
    #Count reads in a FASTQ file (4 lines per read).
    with open(filename) as f:
        return sum(1 for _ in f) // 4

def assembly_stats(contigs_file, min_length=1000):
    count = 0
    total_bp = 0
    for record in SeqIO.parse(contigs_file, "fasta"):
        L = len(record.seq)
        if L > min_length:
            count += 1
            total_bp += L
    return count, total_bp