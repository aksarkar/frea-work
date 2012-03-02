"""Convert regions in an LD neighborhood to a format suitable for plotting

Author: Abhishek Sarkar

Expects BED5 (chr, start, end, annotation, representative) on stdin. Writes CSV
(group, relative start, relative end, annotation) on stdout.
"""
import csv
import itertools as it
import operator as op
import sys

for i, (_, xss) in enumerate(it.groupby(csv.reader(sys.stdin, delimiter='\t'),
                                        key=op.itemgetter(4))):
    for xs in xss:
        print(i, int(xs[1]) - int(xs[4]), int(xs[2]) - int(xs[4]), xs[3], sep=',')
