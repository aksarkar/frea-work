"""Compute expected probability a SNP will be in LD with a feature given the
number of SNPs it is in LD with.

Usage: python exectation.py THRESH

Expects 0-based BED of SNPs with binary annotation in last column on stdin.

Author: Abhishek Sarkar

"""
import csv
import functools
import itertools
import operator
import re
import sys

def intersects(line):
    id_, rest = line.split()
    rest += ';{},1,1'.format(id_)
    inld = (tuple(x.split(',')) for x in rest.split(';'))
    has_r2 = (x for x in inld if len(x) == 3)
    cutoff = (x[0] for x in has_r2 if float(x[2]) >= thresh)
    annots = [snps[x] for x in cutoff if x in snps]
    return len(annots), functools.reduce(operator.or_, annots, 0)

thresh = float(sys.argv[1])
ldtable = sys.argv[2]
snps = {row[3]: int(row[-1]) for row in csv.reader(sys.stdin, delimiter='\t')}
with open(ldtable) as f:
    xs = sorted((intersects(line) for line in f), key=operator.itemgetter(0))
    for k, g in itertools.groupby(xs, key=operator.itemgetter(0)):
        h = [x for _, x in g]
        print(k, sum(h) / len(h))
