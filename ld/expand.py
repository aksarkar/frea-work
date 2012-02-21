"""Expand GWAS SNPs according to LD

Author: Abhishek Sarkar <aksarkar@mit.edu>

Usage: python expand.py THRESH

Expects 0-based BED file on stdin. Writes out 0-based BED4 file with entries as
(chromosome, start, end, pos). The first 3 fields are for the marker in LD; the
last identifies which GWAS SNP they are in LD with.


"""
import csv
import os
import sys

fmt = 'chr{}\t{}\t{}\t{}'

def helper(filename_fmt, index, chrom, pos, thresh, outfile):
    """Look up LD in the specified direction using the appropriate index"""
    if (chrom, pos) not in index:
        return
    with open(filename_fmt.format(chrom)) as f:
        f.seek(index[(chrom, pos)])
        for line in f:
            row = line.split()
            if int(row[0]) != pos:
                return
            if float(row[6]) >= thresh:
                print(fmt.format(chrom, int(row[1]) - 1, row[1], pos))

def lookup(chrom, pos, thresh, outfile):
    """Look up LD in both directions"""
    print(fmt.format(chrom, pos - 1, pos, pos))
    helper(os.path.expanduser('~/hp/hapmap/ld_chr{}_CEU.txt'),
           forward, chrom, pos, thresh, outfile)
    helper(os.path.expanduser('~/hp/hapmap/ld_chr{}_rev.txt'),
           reverse, chrom, pos, thresh, outfile)

if __name__ == '__main__':
    with open(os.path.expanduser('~/hp/hapmap/forward_index.csv')) as f:
        forward = {(int(row[0]), int(row[1])): int(row[2])
                   for row in csv.reader(f)}

    with open(os.path.expanduser('~/hp/hapmap/reverse_index.csv')) as f:
        reverse = {(int(row[0]), int(row[1])): int(row[2])
                   for row in csv.reader(f)}

    thresh = float(sys.argv[1])
    for line in sys.stdin:
        row = line.split()
        chrom = int(row[0][3:])
        pos = int(row[2])
        lookup(chrom, pos, thresh, f)
