"""Permutation test for enrichment of chromatin states

Author: Abhishek Sarkar <aksarkar@mit.edu>

Usage: python3 test_all.py CELLTYPE STATE [NPERMS]

CELLTYPE - cell type
STATE - chromatin state (space-separated list)
NPERMS - number of permutations (default: length of data)

Expects comma-separated annotation file (p-value + state annotations) on stdin.
Writes p-value to stdout.

"""

import csv
import itertools as it
import functools as ft
import multiprocessing as mp
import random
import operator as op
import os
import sys

import bitvector

fmt = '{},"{}",{:.6e},{:.1f}'

def rot(iterable, offset):
    return it.chain(it.islice(iterable, offset, None), 
                    it.islice(iterable, 0, offset))

def perm(ps, annots, cell_type, states, offset):
    return sum(i for i, (p, a) in enumerate(zip(ps, rot(annots, offset))) 
               if a & bitvector.bitvector(states))

def test(pool, data, cell_type, states, nperms):
    ps, annots = zip(*sorted((float(row['p']), int(row[cell_type])) 
                             for row in data))
    f = ft.partial(perm, ps, annots, cell_type, states)
    orig = f(0)
    null = pool.map(f, random.sample(range(len(data)), nperms))
    count = sum(1 if u <= orig else 0 for u in null)
    return fmt.format(cell_type, 
                      ','.join(str(x) for x in states),
                      count / nperms, 1)

if __name__ == '__main__':
    reader = csv.DictReader(sys.stdin)
    data = list(reader)
    pool = mp.Pool(int(os.getenv('LSB_DJOB_NUMPROC')))
    cell_type = sys.argv[1]
    state = [int(x) for x in sys.argv[2].split(' ')]
    nperms = int(sys.argv[3]) if len(sys.argv) > 3 else len(data)
    print(test(pool, data, cell_type, states, nperms))
