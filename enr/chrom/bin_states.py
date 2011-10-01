import csv
import itertools as it
import functools as ft
import math
import sys

import annotate
from bitvector import bitvector as bv

classes = {'promoter': [1, 2, 14],
           'enhancer': [4, 5, 11, 12],
           'insulator': [0],
           'transcribed': [6, 10, 13],
           'repressed': [9],
           'other': [3, 7, 8]}

binsize = int(sys.argv[1])
include = sys.argv[2:]
aggregate = True

r = csv.reader(sys.stdin)
head = next(r)
index_ = [head.index(x) for x in include]
data = [ft.reduce(op.or_, [int(states[i]) for i in index_]) for states in
        sorted(list(r), key=lambda row: float(row[1]))]

expected = [len([x for x in data if x & (1 << state)]) / len(data) 
            for state in range(15)]

count = [0 for _ in range(15)]
w = csv.writer(sys.stdout)
w.writerow(['total', 'feature', 'count', 'expected'])
for i in range(0, len(data), binsize):
    bincounts = (len([x for x in data[i:i + binsize] if x & (1 << state)])
                 for state in range(15))
    count = [c + d for c, d in zip(count, bincounts)]
    if aggregate:
        for k, v in classes.items():
            w.writerow([i + binsize, k, sum(count[j] for j in v), 
                        sum(int((i + binsize) * expected[j]) for j in v)])
    else:
        for j, x in enumerate(count):
            w.writerow([i + binsize, j, x, int((i + binsize) * expected[j])])
