"""Convert regions in an LD neighborhood to a format suitable for plotting

Author: Abhishek Sarkar

Expects BED4 (representative, 0-based start, 0-based end, annotation) on stdin.
Writes (representative, relative start, relative end, annotation) on stdout.
"""
import csv
import itertools as it
import operator as op
import sys

for i, (_, xss) in enumerate(it.groupby(csv.reader(sys.stdin), key=op.itemgetter(0))):
    for xs in xss:
        print(i, int(xs[1]) - int(xs[0]), int(xs[2]) - int(xs[0]), xs[3], sep=',')
