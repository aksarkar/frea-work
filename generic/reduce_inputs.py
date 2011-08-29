"""Generic reducer for enrichment tests

Author: Abhishek Sarkar <aksarkar@mit.edu>

Usage: python reduce_inputs.py

"""

import csv
import sys

def annot(i):
    with open('/broad/shptmp/aksarkar/perm_test.{}'.format(i)) as f:
        r = csv.reader(f)
        w = csv.writer(sys.stdout)
        head = next(r)
        if i == 1:
            w.writerow(['analysis', 'phenotype'] + head)
        w.writerows(lookup[i] + row for row in r)

if __name__ == '__main__':
    with open('/seq/compbio-hp/GWAS/meta/analyses.txt') as f:
        lookup = list(csv.reader(f, delimiter='\t'))
    map(annot, range(1, 4217))
