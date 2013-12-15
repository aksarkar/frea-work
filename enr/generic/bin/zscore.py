import csv
import glob
import itertools
import math
import sys

files = list(itertools.chain.from_iterable(glob.glob(a) for a in sys.argv[1:]))

with open(files[0]) as f:
    data = list(csv.reader(f))
    fold = [float(x[4]) / float(x[5]) for x in data]
    M = fold[:]
    S = [0 for _ in fold]

for i, a in enumerate(files[1:]):
    # Running variance, TAOCP Vol. 2
    with open(a) as f:
        X = [float(x[4]) / float(x[5]) for x in csv.reader(f)]
        N = [m + (x - m) / (i + 2) for m, x in zip(M, X)]
        S = [s + (x - m) * (x - n) for s, x, m, n in zip(S, X, M, N)]
        M = N

for (count, pheno, cell, feature, _, _), f, m, s in zip(data, fold, M, S):
    print(count, pheno, cell, feature, (f - m) / math.sqrt(s), sep=',')
