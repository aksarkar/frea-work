"""Annotate SNPs with enriched pathways

Usage: python annotate.py PATHWAYS

PATHWAYS is expected to be output from DAVID (see
code/enr/david/david.py). Expects one line per SNP on stdin. Writes annotated
matrix on stdout.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import collections
import csv
import signal
import sys

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

enrichments = csv.reader(sys.stdin, delimiter='\t')
edict = {row[3]: set(row[7].split(', ')) for row in enrichments}

header = ['chr', 'start', 'end', 'tag_snp', 'log_p', 'rank', 'snp', 'odds_ratio', 'maf', 'gene_id', 'gene', 'module']
header.extend(k for k in edict)
print(*header, sep='\t')
with open(sys.argv[1]) as f:
    data = (line.split() for line in f)
    for row in data:
        for k in edict:
            row.append(1 if row[9] in edict[k] else 0)
        print(*row, sep='\t')
