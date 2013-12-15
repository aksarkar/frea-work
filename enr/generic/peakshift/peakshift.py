"""Peakshifting permutation test for functional enrichment

Usage: python peakshift.py SNPS REGIONS N THRESH

Shifts regions by random offset in [-THRESH, THRESH] and counts the empirical
frequency over N trials of observing as many overlaps as in the original data.

Expects SNPS and REGIONS to be gzipped, sorted zero-based BED files. Writes
empirical count of trials with as many overlaps, actual count of overlaps, and
mean and variance of overlaps over the shifted data sets on stdout.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import gzip
import itertools
import operator
import random
import sys

def isect(snps, regions):
    region_iter = iter(regions)
    start, end = next(region_iter)
    for snp in snps:
        while end < snp:
            start, end = next(region_iter)
        if snp < start:
            continue
        if end > snp:
            yield snp

with gzip.open(sys.argv[1]) as f:
    snps_raw = (line.split() for line in f)
    snps = {k: sorted(int(x[1]) for x in g) for k, g in
            itertools.groupby(snps_raw, key=operator.itemgetter(0))}
with gzip.open(sys.argv[2]) as f:
    regions_raw = (line.split() for line in f)
    regions = {k: sorted((int(x[1]), int(x[2])) for x in g) for k, g in
               itertools.groupby(regions_raw, key=operator.itemgetter(0)) if k in snps}

ntrials = int(sys.argv[3])
thresh = int(sys.argv[4])
X = sum(len(list(isect(snps[k], regions[k]))) for k in snps if k in regions)
count = 1
running_mean = X
running_variance = 0
offsets = ((-1 if random.random() > .5 else 1) * random.randrange(thresh)
           for _ in itertools.count())
for i in range(ntrials):
    shifted = {k: [(s + o, e + o) for (s, e), o in zip(v, offsets)]
               for k, v in regions.items()}
    Y = sum(len(list(isect(snps[k], shifted[k]))) for k in shifted)
    new_mean = running_mean + (Y - running_mean) / (i + 2)
    running_variance += (Y - running_mean) * (Y - new_mean)
    running_mean = new_mean
    if Y > X:
        count += 1
print('{} {} {:.3f} {:.3f}'.format(count, X, running_mean, running_variance / ntrials))
