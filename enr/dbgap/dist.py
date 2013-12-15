import csv
import itertools
import operator
import sys

types = [int, str, str, str, float]
data = ([f(x) for f, x in zip(types, xs)] for xs in csv.reader(sys.stdin))
for k, g in itertools.groupby(data, key=operator.itemgetter(0, 2, 3)):
    count, cell, annot = k
    xs = sorted(x[4] for x in g)
    print(count, "99th percentile", cell, annot, xs[99 * len(xs) // 100], sep=',')
    print(count, "1st percentile", cell, annot, xs[1 * len(xs) // 100], sep=',')
