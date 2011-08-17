"""Test GWAS p-values for enrichment of chromatin states

Author: Abhishek Sarkar <aksarkar@mit.edu>

Usage: python3 test_enrichment.py [TOP]

TOP - list of top SNPs

Expects comma-separated annotation file (p-value + state annotations)
on stdin. Top SNPs should be listed one rsid per line. Writes p-values
and fold enrichment to stdout.

If TOP is specified, performs Fisher's exact test on top
SNPs.

"""
import csv
import itertools as it
import random
import operator as op
import sys

import scipy.stats

import bitvector

fmt = '{},"{}",{:.6e},{:.1f}'

def test_top(data, cell_type, states, top):
    xs = [row['rsid'] for row in data if 
          int(row[cell_type]) & bitvector.bitvector(states)]
    xa = len([x for x in xs if x in top])
    xb = len(xs) - xa
    ya = len(top) - xa
    yb = len(data) - len(xs) - ya
    obs = xa / (xa + ya)
    exp = len(xs) / len(data)
    return fmt.format(cell_type, 
                      ','.join(str(x) for x in states),
                      scipy.stats.fisher_exact([[xa, xb], [ya, yb]])[1],
                      obs / exp)
def rot(iterable, offset):
    return it.chain(it.islice(iterable, offset, None), 
                    it.islice(iterable, 0, offset))

def perm_helper(ps, annots, cell_type, states):
    return sum(i for i, (p, a) in enumerate(zip(ps, annots)) 
               if a & bitvector.bitvector(states))

def perm(data, cell_type, states, nperms=None):
    if nperms is None:
        nperms = len(data)
    ps, annots = zip(*sorted((float(row['p']), int(row[cell_type])) for row in data))
    orig = perm_helper(ps, annots, cell_type, states)
    count = sum(1 if perm_helper(ps, rot(annots, i), cell_type, states) <= orig 
                else 0 for i in random.sample(range(len(data)), nperms))
    return fmt.format(cell_type, 
                      ','.join(str(x) for x in states),
                      count / nperms, 1)

def test(data, cell_type, states, top=[]):
    """Test for enrichment of p-values for given cell_type and states"""
    if top:
        return test_top(data, cell_type, states, top)
    else:
        return perm(data, cell_type, states, 1000)

if __name__ == '__main__':
    reader = csv.DictReader(sys.stdin)
    data = list(reader)
    top = []
    if len(sys.argv) > 1:
        with open(sys.argv[1]) as f:
            top = set([x.strip() for x in f])
    print('cell_type,state,p,fold')
    candidates = [[11, 12], [4, 5], [11], [12], [4], [5]]
    print('\n'.join(test(data, cell_type, states, top)
                    for cell_type in reader.fieldnames[2:3]
                    for states in candidates[:1]))
