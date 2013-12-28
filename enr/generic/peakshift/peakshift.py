"""Peakshifting permutation test for functional enrichment

Usage: python peakshift.py SNPS REGIONS N THRESH

Shifts regions by random offset in [-THRESH, THRESH] and counts the empirical
frequency over N trials of observing as many overlaps as in the original data.

Expects SNPS and REGIONS to be sorted zero-based BED files. Writes
empirical count of trials with as many overlaps, actual count of overlaps, and
mean and variance of overlaps over the shifted data sets on stdout.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import bz2
import gzip
import itertools
import operator
import os
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

def handle(filename, loader):
    def helper(fn):
        with fn(filename, 'rt') as f:
            return loader(f)
    if filename == '-':
        return loader(sys.stdin)
    _, ext = os.path.splitext(filename)
    if ext == '.bz2':
        return helper(bz2.open)
    elif ext == '.gz':
        return helper(gzip.open)
    else:
        return helper(open)

def load_snps(f):
    snps_raw = (line.split() for line in f)
    return {k: [int(x[1]) for x in g] for k, g in
            itertools.groupby(snps_raw, key=operator.itemgetter(0))}

def load_regions(f):
    regions_raw = (line.split() for line in f)
    return {k: [(int(x[1]), int(x[2])) for x in g] for k, g in
            itertools.groupby(regions_raw, key=operator.itemgetter(0)) if k in snps}

def trial(snps, regions, offsets):
    shifted = {k: [(s + o, e + o) for (s, e), o in zip(v, offsets)]
               for k, v in regions.items()}
    return sum(len(list(isect(snps[k], shifted[k]))) for k in shifted)

if sys.argv[1] == sys.argv[2] == '-':
    print('error: SNPS and REGIONS cannot both be -')
    sys.exit(1)

snps = handle(sys.argv[1], load_snps)
regions = handle(sys.argv[2], load_regions)

ntrials = int(sys.argv[3])
if ntrials <= 0:
    print('error: N must be positive')
    sys.exit(1)
thresh = int(sys.argv[4])
if thresh <= 0:
    print('error: thresh must be positive')
    sys.exit(1)
X = sum(len(list(isect(snps[k], regions[k]))) for k in snps if k in regions)
offsets = ((-1 if random.random() > .5 else 1) * random.randrange(thresh)
           for _ in itertools.count())
Y = trial(snps, regions, offsets)
count = 1
if Y > X:
    count += 1
running_mean = Y
running_variance = 0
for i in range(ntrials - 1):
    Y = trial(snps, regions, offsets)
    new_mean = running_mean + (Y - running_mean) / (i + 2)
    running_variance += (Y - running_mean) * (Y - new_mean)
    running_mean = new_mean
    if Y > X:
        count += 1
print('{} {} {:.3f} {:.3f}'.format(count, X, running_mean, running_variance))
