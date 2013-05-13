"""Filter WTCCC1 genotypes according to lifted over Affy probes

Usage: python genotypes.py PROBES GENOTYPES

Expects input in IMPUTE2 format, outputs in IMPUTE2 format using positions,
rsids in PROBES.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import io
import gzip
import sys

with open(sys.argv[1]) as f:
    rows = (line.split() for line in f)
    next(rows)
    snps = {row[3]: row for row in rows}
rows = (line.split() for line in sys.stdin)
b137 = (snps[row[1]] + row for row in rows if row[1] in snps)
for row in b137:
    print(row[6], row[3], row[1] + 1, row[5], row[6], *row[11:])
