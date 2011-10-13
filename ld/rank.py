"""Count how many independent SNPs there are.

Author: Abhishek Sarkar <aksarkar@mit.edu>

Usage: python rank.py RANKLIST ANNOT

RANKLIST - BED file of GWAS markers, ascending p-value
ANNOT - BED file of SNPs in LD with GWAS markers
"""
import csv
import itertools as it
import operator as op
import sys

with open(sys.argv[1]) as f:
    ranklist = set((row[0], row[1]) for row in csv.reader(f, delimiter='\t'))

with open(sys.argv[2]) as f:
    count = 0
    for k, g in it.groupby(csv.reader(f, delimiter='\t'), key=op.itemgetter(3)):
        if not ranklist:
            break
        else:
            t = set((row[0], row[1]) for row in g)
            if ranklist & t:
                count += 1
            ranklist -= t
    count += len(ranklist)
    print(count)
