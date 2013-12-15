import csv
import itertools
import operator
import sys

types = [int, str, str, str, float, float]
data = ([f(x) for f, x in zip(types, row)] for row in csv.reader(sys.stdin))
for k, g in itertools.groupby(data, key=operator.itemgetter(1)):
    x = list(g)
    c = x[-1][4] + 1
    dev = (bin_[:4] + [(bin_[4] - bin_[5]) / c] for bin_ in x if bin_[0] <= 150000)
    for r in dev:
        print(*r, sep=',')
