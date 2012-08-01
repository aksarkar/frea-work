"""Test GWAS p-values for enrichment of DNAse across cell types

Author: Abhishek Sarkar <aksarkar@mit.edu>

Usage: python3 test_dnase_enrichment.py NSNPS CELLTYPES [TOP]

NSNPS - total number of SNPs
CELLTYPES - list of cell types to test for enrichment
TOP - list of top SNPs

Expects annotation file containing (rsid, p value, cell type) triples
on stdin. Cell types and top SNPs should be listed one per line. Writes
p-values and fold enrichment to stdout.

"""
import csv
import itertools as it
import operator as op
import sys

import scipy.stats

fmt = '{},{:.6e},{:.1f}'

def test(cell_type, data, top, nsnps):
    xs = [x for x in data if x[2] == cell_type]
    if not xs:
        return fmt.format(cell_type, 1, 0)
    elif top:
        xa = len([x for x in xs if x[0] in top])
        xb = len(xs) - xa
        ya = len(top) - xa
        yb = nsnps - len(xs) - ya
        obs = xa / len(top)
        exp = len(xs) / len(data)
        return fmt.format(cell_type, 
                          scipy.stats.fisher_exact([[xa, xb], [ya, yb]])[1],
                          obs / exp)
    else:
        raise NotImplementedError('RR test')

if __name__ == '__main__':
    r = csv.reader(sys.stdin)
    head = next(r)
    types = [str, float, str]
    data = [[f(x) for f, x in zip(types, row)] for row in r]
    nsnps = int(sys.argv[1])
    with open(sys.argv[2]) as f:
        cell_types = [x.strip() for x in f]
    with open(sys.argv[3]) as f:
        top = set([x.strip() for x in f])
    print('cell_type,p,fold')
    print('\n'.join(test(cell_type, data, top, nsnps)
                    for cell_type in cell_types))
