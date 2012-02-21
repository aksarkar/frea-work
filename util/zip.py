"""Concatenate annotations in a format suitable for visualization

Author: Abhishek Sarkar <aksarkar@mit.edu>

Usage: python zip.py F1 [F2 ...]

Expects (group, start, end, annotation) in each input file. Writes (y, start,
end, annotation) on stdout.
"""
import csv
import itertools as it
import operator as op
import sys

gss = [it.groupby(csv.reader(open(f)), key=op.itemgetter(0))
       for f in sys.argv[1:]]
i = 0
for gs in zip(*gss):
    for _, g in gs:
        for xs in g:
            print(i, *xs[1:], sep=',')
        i += 1
