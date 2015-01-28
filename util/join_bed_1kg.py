import gzip
import itertools
import operator
import sys

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

def output(bed_entry, ref_entry):
    # chr_, _, _, _, score = bed_entry
    chr_, _, b0, b1, score = bed_entry
    name, pos, a0, a1 = ref_entry
    if len(a0) > 1:
        a0 = 'I{}'.format(len(a0))
        a1 = 'D'
    elif len(a1) > 1:
        a0 = 'D'
        a1 = 'I{}'.format(len(a1))
    if (b0, b1) != (a0, a1) and (b1, b0) != (a0, a1):
        print(chr_, pos, b0, b1, a0, a1, file=sys.stderr)
        return
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
    print(chr_, pos - 1, pos + len(a0) - 1, '|'.join([name, chr_, str(end), str(end), delta]), score, sep='\t')
    
if __name__ == '__main__':
    import math
    rawdata = (line.split() for line in sys.stdin)
    next(rawdata)
    # bed = parse([str, int, int, str, float], (line.split() for line in sys.stdin))
    bed = ([row[0], int(row[4]), row[2], row[3], -math.log(float(row[8]), 10)] for row in rawdata)
    for k, g in itertools.groupby(bed, key=operator.itemgetter(0)):
        with gzip.open('/broad/compbio/aksarkar/data/1kg/ALL.{}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.legend.gz'.format(k), 'rt') as f:
            next(f)
            ref = parse([str, int, str, str], (line.split() for line in f))
            for pair in join(g, ref, operator.itemgetter(1)):
                output(*pair)
