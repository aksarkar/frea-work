"""Filter WTCCC1 genotypes according to lifted over Affy probes

Usage: python genotypes.py PROBES OUTDIR [CHR1 ...]

Expects input in IMPUTE2? format, outputs in IMPUTE2 format using positions,
rsids in PROBES.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import gzip
import os
import sys

with open(sys.argv[1]) as f:
    rows = (line.split() for line in f)
    probes = {row[3]: row for row in rows}

basedir = sys.argv[2]

for chromosome in sys.argv[3:]:
    with gzip.open(chromosome) as f, open(chromosome[:-3], 'w') as g:
        decode = (str(line, 'utf-8') for line in f)
        rows = (line.split() for line in decode)
        b137 = (probes[row[0]] + row for row in rows if row[0] in probes)
        for row in b137:
            print(row[3], row[4], row[1], row[10:], file=g)
