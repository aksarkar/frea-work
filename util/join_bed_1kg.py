""" Join SNP lists to 1KG reference information

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import collections
import gzip
import itertools
import operator
import math
import sys

bedtools_format = [str, int, int, str, float]
plink_format = [lambda x: 'chr{}'.format(x), str, int, int, str, str]
legend_format = [str, int, str, str, str]

def parse(types, iterable):
    for row in iterable:
        yield [f(x) for f, x in zip(types, row)]

def join(bed, ref, key, key2=None):
    if key2 is None:
        key2 = key
    bed_buckets = itertools.groupby(bed, key)
    ref_buckets = itertools.groupby(ref, key2)
    k1, g1 = next(bed_buckets)
    k2, g2 = next(ref_buckets)
    while True:
        if k1 == k2:
            for pair in itertools.product(g1, g2):
                yield pair
            k1, g1 = next(bed_buckets)
            k2, g2 = next(ref_buckets)
        elif k1 < k2:
            k1, g1 = next(bed_buckets)
        else:
            k2, g2 = next(ref_buckets)

def get_pouyak_name(chr_, name, pos, a0, a1):
    if len(a0) <= len(a1):
        start = pos
        end = pos
        if len(a0) == len(a1):
            delta = a1
        else:  # insertion
            delta = a1[len(a0):]
    else:  # deletion
        start = pos + len(a1)
        end = pos + len(a0) - len(a1)
        delta = ""
    return '|'.join([name, chr_, str(end), str(end), delta])

def get_plink_alleles(a0, a1):
    if len(a0) > 1:
        a0 = 'I{}'.format(len(a0))
        a1 = 'D'
    elif len(a1) > 1:
        a0 = 'D'
        a1 = 'I{}'.format(len(a1))
    return a0, a1

def is_alignable(a0, a1, b0, b1):
    return (b0, b1) == (a0, a1) or (b1, b0) == (a0, a1)


def lookup(input_file, input_format, output_fn, input_sort_key=operator.itemgetter(0),
         input_join_key=operator.itemgetter(1), ref_join_key=operator.itemgetter(1)):
    """Lookup input SNPs in 1KG

    input_file - file-like object
    input_format - list of column types
    output_fn - function to output hits
    input_sort_key - index of chromosome column
    input_join_key - key to lookup in 1KG
    ref_join_key - key to match in 1KG
    """
    data = (line.split('\t') for line in input_file)
    bed = parse(input_format, data)
    for k, g in itertools.groupby(bed, key=input_sort_key):
        with gzip.open('/broad/compbio/aksarkar/data/1kg/ALL.{}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.legend.gz'.format(k), 'rt') as f:
            next(f)
            ref = parse(legend_format, (line.split() for line in f))
            keep = (x for x in ref if x[-1] == 'SNP')
            for pair in join(g, keep, input_join_key, ref_join_key):
                output_fn(*pair)

if __name__ == '__main__':
    def output(bed_entry, ref_entry, check_strand=False):
        _, chr_, _, score = bed_entry
        name, pos, a0, a1, _ = ref_entry
        if score == 0:
            print(chr_, pos, file=sys.stderr)
            return
        if check_strand and not is_alignable(a0, a1, b0, b1):
            print(chr_, pos, b0, b1, a0, a1, file=sys.stderr)
            return
        print(chr_, pos - 1, pos + len(a0) - 1, get_pouyak_name(chr_, name, pos, a0, a1), -math.log10(score), sep='\t')
    with gzip.open('/broad/compbio/aksarkar/data/gwas-summary-stats/t1d_bradfield_10.1371_journal.pgen.1002293/hg19_gwas_t1d_bradfield_4_18_0.gz', 'rt') as f:
        next(f)
        lookup(f, [str, lambda x: 'chr{}'.format(x), int, float], output, input_sort_key=operator.itemgetter(1), input_join_key=operator.itemgetter(2))
