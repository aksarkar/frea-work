"""Permute elements between cell types

Usage: permute-regions.py NOUTPUTS

Expects list of BED files on stdin (all must have the same number of fields).

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import csv
import gzip
import itertools
import math
import random
import sys

def load(filename):
    with gzip.open(filename) as f:
        return [str(line, 'utf-8').strip() for line in f]

noutputs = int(sys.argv[1])
data = [load(filename.strip()) for filename in sys.stdin]
dist = [len(x) for x in data]
elems = list(itertools.chain.from_iterable(data))
for i in range(noutputs):
    with open('{:0{w}d}.bed'.format(i, w=int(math.log10(noutputs))),
              'w') as f:
        for j in random.sample(elems, random.choice(dist)):
            print(j, file=f)
