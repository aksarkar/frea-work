"""Generate distribution of functional elements around tagged markers

Author: Abhishek Sarkar

Expects (group, start, end, annotation) tuples on stdin. Writes (bin,
annotation, count) on stdout.

"""
import collections
import csv
import itertools as it
import operator as op
import sys

def bins(start, end, binsize, annot):
    return [(x, annot) for x in range((start // binsize) * binsize,
                                      (end // binsize) * binsize + binsize,
                                      binsize)]

binsize = int(sys.argv[1])
types = [int, int, int, str]
data = ([f(x) for f, x in zip(types, xs)] for xs in csv.reader(sys.stdin))
result = collections.Counter()
for _, xss in it.groupby(data, key=op.itemgetter(0)):
    result.update(it.chain.from_iterable(bins(start, end, binsize, annot)
                                         for (_, start, end, annot) in xss))
for (bin, annot), count in sorted(result.items()):
    print(bin, annot, count, sep=',')
