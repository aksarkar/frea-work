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
import collections
import gzip
import itertools
import operator
import os
import random
import sys

import scipy.stats

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
    head = collections.defaultdict(list)
    tail = collections.defaultdict(list)
    for k, g in itertools.groupby(snps_raw, key=operator.itemgetter(0)):
        parsed = ((int(pos), float(score)) for _, pos, _, _, score in g)
        for p, s in parsed:
            if s >= 2.214570:
                head[k].append(p)
            else:
                tail[k].append(p)
    return head, tail

def load_regions(f):
    regions_raw = (line.split() for line in f)
    return {k: [(int(x[1]), int(x[2])) for x in g] for k, g in
            itertools.groupby(regions_raw, key=operator.itemgetter(0))}

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

def noverlaps(snps, regions):
    inside, outside = 0, 0
    for k in snps:
        inside += len(list(isect(snps[k], regions.get(k, []))))
        outside += len(snps[k]) - inside
    return inside, outside

def trial(snps, regions, offsets=None):
    if offsets is not None:
        regions = {k: [(s + o, e + o) for (s, e), o in zip(v, offsets)]
                   for k, v in regions.items()}
    head, tail = snps
    return noverlaps(head, regions)[0], noverlaps(tail, regions)[0]

def test(snps, regions, thresh):
    X = trial(snps, regions)
    offsets = ((-1 if random.random() > .5 else 1) * random.randrange(thresh)
               for _ in itertools.count())
    running_mean = trial(snps, regions, offsets)
    for i in range(ntrials - 1):
        Y = trial(snps, regions, offsets)
        running_mean = [r + (y - r) / (i + 2) for r, y in zip(running_mean, Y)]
    return scipy.stats.fisher_exact([X, running_mean], alternative='greater')

if __name__ == '__main__':
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
    print('{:.3f} {:.3f}'.format(*test(snps, regions, thresh)))
