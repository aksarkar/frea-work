"""Take the union of annotations over LD blocks

Author: Abhishek Sarkar <aksarkar@mit.edu>

Usage: python group.py FN

FN - the function used to combine annotations. It must be one of 'state' or
'union'.

Reads BED file on stdin. The fourth field identifies the GWAS marker the
current SNP is in LD with. The last field is the annotation for the current
SNP.

"""
import csv
import itertools as it
import functools as ft
import operator as op
import sys

def state(xs):
    return ft.reduce(op.or_, (1 << int(x[-1]) for x in xs))

def union(xs):
    return ft.reduce(op.or_, (int(x[-1]) for x in xs))

fn = {'state': state, 'union': union}[sys.argv[1]]

for k, g in it.groupby(csv.reader(sys.stdin, delimiter='\t'),
                       key=op.itemgetter(3)):
    xs = list(g)
    print('\t'.join(xs[0][:4]), fn(xs), sep='\t')
