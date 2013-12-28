"""Compute hypergeometric p-values for functional enrichment with MAF matched
negative sets

Usage: python fisher.py MAFS THRESH NTRIALS

Expects BED file of SNPs with score = negative log-transformed p-values on
stdin. MAF is space-separated (rsid, MAF) pairs.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import bisect
import collections
import gzip
import itertools
import math
import random
import sys

import scipy.stats

thresh = float(sys.argv[2])
ntrials = int(sys.argv[3])
breaks = [0.0005, 0.001, 0.005, 0.01, 0.03, 0.05, 0.08, 0.10, 0.15, 0.20,
          0.30, 0.40, 0.50]

snps_raw = (line.split() for line in sys.stdin)
snps = {row[3]: (float(row[4]), int(row[5])) for row in snps_raw}

with gzip.open(sys.argv[1]) as f:
    data = (str(l, 'utf-8').split() for l in f)
    filter_ = (d for d in data if d[0] in snps)
    overlap_bins = collections.defaultdict(list)
    nonoverlap_bins = collections.defaultdict(list)
    for k, g in itertools.groupby(filter_, key=lambda x: bisect.bisect(breaks, float(x[1]))):
        for rsid, _ in g:
            if snps[rsid][1]:
                overlap_bins[k].append(rsid)
            else:
                nonoverlap_bins[k].append(rsid)
overlaps = {k: len(overlap_bins[k]) for k in overlap_bins}
total_overlaps = sum(overlaps.values())

I = itertools.chain.from_iterable
num_top_overlaps = 0
for s in snps:
    if snps[s][0] > thresh and snps[s][1]:
        num_top_overlaps += 1
for k in overlaps:
    assert k in nonoverlap_bins
    assert nonoverlap_bins[k]
num_top_matches = 0
for _ in range(ntrials):
    matches = I((random.choice(nonoverlap_bins[k]) for _ in range(v)) for k, v in overlaps.items())
    for s in matches:
        if snps[s][0] > thresh:
            num_top_matches += 1
num_top_matches /= ntrials
print(*scipy.stats.fisher_exact([[num_top_overlaps, total_overlaps - num_top_overlaps],
                                 [num_top_matches, total_overlaps - num_top_matches]],
                                 alternative='greater'))
