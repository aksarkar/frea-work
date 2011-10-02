"""Generic reducer for enrichment tests

Author: Abhishek Sarkar <aksarkar@mit.edu>

Usage: python reduce_inputs.py

"""

import csv
import os
import sys

def annot(i):
    with open('/broad/shptmp/aksarkar/map_inputs.{}'.format(i + 1)) as f:
        r = csv.reader(f)
        w = csv.writer(sys.stdout)
        w.writerows(lookup[i] + row for row in r)

if __name__ == '__main__':
    with open('/seq/compbio-hp/GWAS/meta/analyses.txt') as f:
        lookup = list(csv.reader(f, delimiter='\t'))
    for i in range(4216):
        annot(i)