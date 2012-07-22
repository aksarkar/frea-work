"""Expand GWAS SNPs according to Hapmap 3 LD

Author: Abhishek Sarkar <aksarkar@mit.edu>

Usage: python expand.py THRESH

Expects 0-based BED file with hg18 coordinates on stdin. Writes out 0-based
BED5 file (chromosome, start, end, representative, p, R²).

"""
import csv
import itertools
import functools
import operator
import os
import sys

def output(chrom, pos, rep, p, corr):
    print(chrom, pos - 1, pos, rep, p, corr, sep='\t')

def helper(markers, ldfile, index):
    """Look up LD in the specified direction using the appropriate index"""
    with open(ldfile) as f:
        for chrom, _, pos, _, p in markers:
            output(chrom, pos, pos, p, 1)
            if (int(chrom[3:]), pos) in index:
                f.seek(index[(int(chrom[3:]), pos)])
                # pos1, pos2, pop, rsid1, rsid2, D', R², LOD
                types = [int, int, str, str, str, float, float, float]
                rows = ([g(x) for g, x in zip(types, line.split())] for line in f)
                inld = itertools.takewhile(lambda row: row[0] == pos, rows)
                meets_thresh = filter(lambda row: row[6] >= thresh, inld)
                for row in meets_thresh:
                    output(chrom, row[1], pos, p, row[6])

def lookup(chrom, markers, findex, rindex):
    """Look up LD in both directions"""
    h = functools.partial(helper, list(markers))
    h('/broad/compbio/aksarkar/hapmap/ld_{}_CEU.txt'.format(chrom), findex)
    h('/broad/compbio/aksarkar/hapmap/ld_{}_rev.txt'.format(chrom), rindex)

if __name__ == '__main__':
    with open('/broad/compbio/aksarkar/hapmap/forward_index.csv') as f:
        findex = {(int(row[0]), int(row[1])): int(row[2])
                  for row in csv.reader(f)}
    with open('/broad/compbio/aksarkar/hapmap/reverse_index.csv') as f:
        rindex = {(int(row[0]), int(row[1])): int(row[2])
                  for row in csv.reader(f)}
    thresh = float(sys.argv[1])
    key = operator.itemgetter(0)
    types = [str, int, int, str, float]
    markers = ([g(x) for g, x in zip(types, row)]
               for row in csv.reader(sys.stdin, delimiter='\t'))
    for k, g in itertools.groupby(sorted(markers, key=key), key=key):
        lookup(k, g, findex, rindex)
