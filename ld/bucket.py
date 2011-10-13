"""Bucket GWAS SNPs which are in LD with each other

Author: Abhishek Sarkar <aksarkar@mit.edu>

Usage: python rank.py GWAS EXPANDED

Writes BED file to stdout. The fourth field identifies the GWAS marker which is
the representative for the current marker.

"""

import csv
import itertools as it
import operator as op
import sys

with open(sys.argv[1]) as f:
    gwas = {(row[0], int(row[1])): (row[0], int(row[1]), float(row[4]))
            for row in csv.reader(f, delimiter='\t')}

with open(sys.argv[2]) as f:
    for k, g in it.groupby(csv.reader(f, delimiter='\t'), key=op.itemgetter(3)):
        inld = [(row[0], int(row[1])) for row in g]
        bucket = [b for b in inld if b in gwas]
        (_, _, p), (c, pos) = min((gwas[b], b) for b in bucket)
        for b in bucket:
            gwas[b] = (c, pos, p)

for c, p in sorted(gwas.keys()):
    print(c, p, p + 1, '{},{}\t{}'.format(*gwas[(c, p)]), sep='\t')
