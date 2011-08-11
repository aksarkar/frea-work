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
import sys

import scipy.stats

import annotate
import bitvector

fmt = '{},"{}",{:.6e},{:.1f}'

def test(data, cell_type, states, top=[]):
    """Test for enrichment of p-values for given cell_type and states"""
    mask = bitvector.bitvector(states)
    xs = [row for row in data if int(row[cell_type]) & mask]
    ys = [row for row in data if not int(row[cell_type]) & mask]
    if not xs or not ys:
        return fmt.format(cell_type, ','.join(str(x) for x in states), 1, 0)
    elif top:
        xa = len([x for x in xs if x['rsid'] in top])
        xb = len(xs) - xa
        ya = len([y for y in ys if y['rsid'] in top])
        yb = len(ys) - ya
        obs = xa / (xa + ya)
        exp = len(xs) / len(data)
        return fmt.format(cell_type, 
                          ','.join(str(x) for x in states),
                          scipy.stats.fisher_exact([[xa, xb], [ya, yb]])[1],
                          obs / exp)
    else:
        raise NotImplementedError('RR test')

if __name__ == '__main__':
    data = list(csv.DictReader(sys.stdin))
    top = []
    if len(sys.argv) > 1:
        with open(sys.argv[1]) as f:
            top = set([x.strip() for x in f])
    print('cell_type,state,p,fold')
    candidates = [[11, 12], [4, 5], [11], [12], [4], [5]]
    print('\n'.join(test(data, cell_type, states, top)
                    for cell_type in annotate.celltypes
                    for states in candidates))
