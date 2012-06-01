"""Create a matrix of genotypes for WTCCC data

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import codecs
import collections
import csv
import gzip
import functools
import itertools
import os
import sys

def wtccc_convert(row):
    work = [float(x) for x in row if x]
    r = round
    return [int(0 * r(x) + r(y) + 2 * r(z)) for x, y, z in
            zip(work[::3], work[1::3], work[2::3])]

def wtccc_load_chromosome(hits, infile, outfile):
    decode = codecs.iterdecode(gzip.open(infile), 'ascii')
    filter_null = (line for line in decode if '\x00' not in line)
    for row in csv.reader(filter_null, delimiter=' '):
        if len(row) > 5 and row[1] in hits:
            print(row[1], *wtccc_convert(row[5:]), file=outfile)

def t1dgc_load_snp(exclude, rsid, infile, outfile):
    with gzip.open(infile) as f:
        decode = codecs.iterdecode(f, 'ascii')
        print(rsid, *[int(row[2]) - 1 for row in csv.reader(decode, delimiter='\t')
                      if row[1] not in exclude],
              file=outfile)

def t1dgc_load_chromosome(exclude, hits, basedir, outfile):
    for k in hits:
        snpfile = os.path.join(basedir, 'control-{}.gz'.format(k))
        if os.path.exists(snpfile):
            t1dgc_load_snp(exclude, k, snpfile, outfile)
        else:
            print('Could not find {}'.format(snpfile))

def load(hits, load_chromosome, pattern, outfile):
    for k in hits:
        load_chromosome(hits[k], pattern.format(k), outfile)

if __name__ == '__main__':
    with gzip.open('/broad/compbio/aksarkar/t1d/data/30k/random.bed.gz') as f:
        hits = collections.defaultdict(set)
        for row in csv.reader(codecs.iterdecode(f, 'ascii'), delimiter='\t'):
            hits[int(row[0][3:])].add(row[3])

    L = functools.partial(load, hits, wtccc_load_chromosome)
    with open('/broad/compbio/aksarkar/t1d/enhancers/genotypes/random/cases', 'w') as f:
        L('/broad/compbio/lward/GWAS/ebi_ega_downloads/WTCCC_7_Diseases/'
          'T1D/imputed/T1D_chr{:02d}_imputed.txt.gz', f)

    with open('/broad/compbio/aksarkar/t1d/enhancers/genotypes/random/controls_1958bc', 'w') as f:
        L('/broad/compbio/lward/GWAS/ebi_ega_downloads/1958BC/'
          'WTCCC1_1958BC/imputed/58C_chr{:02d}_imputed.txt.gz', f)

    with open('/broad/compbio/aksarkar/t1d/enhancers/genotypes/random/controls_nbs', 'w') as f:
        L('/broad/compbio/lward/GWAS/ebi_ega_downloads/NBS/'
          'WTCCC1/imputed/NBS_chr{:02d}_imputed.txt.gz', f)

    with open('/broad/compbio/aksarkar/wtccc/EGAD00000000030/'
              'T1DGC-GWAS-2008-EGA-Sample-Exclusions-2009-05-08.txt') as f:
        exclude = set(line.strip() for line in f)
        
    M = functools.partial(load, hits, functools.partial(t1dgc_load_chromosome, exclude))
    with open('/broad/compbio/aksarkar/t1d/enhancers/genotypes/random/controls_t1dgc', 'w') as f:
        M('/broad/compbio/aksarkar/wtccc/EGAD00000000030/{:02d}', f)
