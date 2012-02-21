"""Generate LD neighborhoods for visualization

Author: Abhishek Sarkar <aksarkar@mit.edu>

Expects LD-expanded BED4 (chr, 0-based pos, pos+1, representative) on stdin.
Writes (chr, start, end, representative) on stdout.
"""
import csv
import itertools as it
import operator as op
import sys

types = [str, int, int, int]
data = [[f(x) for f, x in zip(types, xs)] for xs in
        csv.reader(sys.stdin, delimiter='\t')]

for _, xs in it.groupby(data, op.itemgetter(3)):
    ys = list(xs)
    lo = min(linked for _, linked, _, pos in ys)
    hi = max(linked for _, linked, _, pos in ys)
    chr_, _, _, pos = ys[0]
    print(chr_, lo, hi, pos, sep='\t')
