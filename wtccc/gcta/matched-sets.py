import collections
import gzip
import itertools
import random
import sys

def build_bins(snps, table):
    overlap_bins = collections.defaultdict(list)
    nonoverlap_bins = collections.defaultdict(list)
    for rsid, *key in table:
        k = tuple(key)
        if rsid in snps:
            overlap_bins[k].append(rsid)
        else:
            nonoverlap_bins[k].append(rsid)
    return overlap_bins, nonoverlap_bins

def match(overlap_counts, nonoverlap_bins):
    I = itertools.chain.from_iterable
    return I((random.choice(nonoverlap_bins[k]) for _ in range(v)) for k, v in
             overlap_counts.items() if k in nonoverlap_bins)

if __name__ == '__main__':
    with open(sys.argv[1]) as f:
        snps = set(line.strip() for line in f)

    with gzip.open(sys.argv[2]) as f:
        data = (str(l, 'utf-8').split() for l in f)
        overlap_bins, nonoverlap_bins = build_bins(snps, data)

    random.seed(int(sys.argv[3]))

    overlap_counts = {k: len(overlap_bins[k]) for k in overlap_bins}
    for s in match(overlap_counts, nonoverlap_bins):
        print(s)
