"""Generate matrix of WTCCC genotypes for SNPs on a single chromosome

Usage: python3 load.py CHROMOSOME

Expects list of rsids on stdin (one per line)

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import codecs
import collections
import csv
import gzip
import functools
import os
import sys

def argmax(iterable):
    return max((x, i) for i, x in enumerate(iterable))[1]

def snptest_call(row):
    work = [float(x) for x in row if x]
    return [argmax(xs) for xs in
            zip(work[::3], work[1::3], work[2::3])]

def load_chromosome(hits, infile, outfile):
    with gzip.open(infile) as f:
        decode = codecs.iterdecode(f, 'ascii')
        filter_null = (line for line in decode if '\x00' not in line)
        for row in csv.reader(filter_null, delimiter=' '):
            if len(row) > 5 and row[1] in hits:
                print(row[1], *snptest_call(row[5:]), file=outfile)

if __name__ == '__main__':
    i = int(sys.argv[1])
    hits = set(line.strip() for line in sys.stdin)
    if not hits:
        sys.exit(0)
    basedir = os.getcwd()
    J = os.path.join
    L = functools.partial(load_chromosome, hits)
    with open(J(basedir, 'cases.{:02d}'.format(i)), 'w') as f:
        L('/broad/compbio/lward/GWAS/ebi_ega_downloads/WTCCC_7_Diseases/'
          'T1D/imputed/T1D_chr{:02d}_imputed.txt.gz'.format(i), f)

    with open(J(basedir, 'controls_1958bc.{:02d}'.format(i)), 'w') as f:
        L('/broad/compbio/lward/GWAS/ebi_ega_downloads/1958BC/'
          'WTCCC1_1958BC/imputed/58C_chr{:02d}_imputed.txt.gz'.format(i), f)

    with open(J(basedir, 'controls_nbs.{:02d}'.format(i)), 'w') as f:
        L('/broad/compbio/lward/GWAS/ebi_ega_downloads/NBS/'
          'WTCCC1/imputed/NBS_chr{:02d}_imputed.txt.gz'.format(i), f)
