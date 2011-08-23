"""Fisher's exact test for enrichment of chromatin states

Author: Abhishek Sarkar <aksarkar@mit.edu>

Usage: python3 test_enrichment.py [TOP]

TOP - list of top SNPs

Expects comma-separated annotation file (p-value + state annotations)
on stdin. Top SNPs should be listed one rsid per line. Writes p-values
and fold enrichment to stdout.

"""
import csv
import itertools as it
import sys

import scipy.stats

import bitvector

fmt = '{},"{}",{:.6e},{:.1f}'

def test(data, cell_type, states, top):
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

if __name__ == '__main__':
    reader = csv.DictReader(sys.stdin)
    data = list(reader)
    with open(sys.argv[1]) as f:
        top = set([x.strip() for x in f])
    print('cell_type,state,p,fold')
    candidates = [[4, 5, 11, 12], [11, 12], [4, 5], [11], [12], [4], [5]]
    print('\n'.join(test(data, cell_type, states, top)
                     for cell_type in reader.fieldnames[2:]
                     for states in candidates))
    
