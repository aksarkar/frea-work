"""Fix representatives when lifting LD-expanded BED files from hg18 to hg19

Author: Abhishek Sarkar <aksarkar@mit.edu>

Expects LD-expanded BED4 on stdin. Writes BED4 on stdout.
"""
import csv
import itertools as it
import operator as op
import sys

for k, xss in it.groupby(csv.reader(sys.stdin, delimiter='\t'), key=op.itemgetter(3)):
    yss = list(xss)
    for ys in yss:
        ys[3] = yss[0][1]
        print(*ys, sep='\t')
