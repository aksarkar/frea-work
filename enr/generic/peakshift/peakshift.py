"""Peakshifting permutation test for functional enrichment

Usage: python peakshift.py SNPS REGIONS N THRESH MINSCORE

Shift regions by random offset in [-THRESH, THRESH] and count the empirical
frequency over N trials of observing as many loci with score greater than
MINSCORE overlapping as in the original data.

Expects SNPS and REGIONS to be sorted zero-based BED files. Writes to stdout:

- observed count of overlaps
- mean of count of overlaps over the shifted data sets
- variance ...
- Kolmogorov--Smirnov p-value for normality
- empirical p-value

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import bz2
import collections
import gzip
import itertools
import operator
import os
import math
import random
import signal
import sys

import scipy.stats

def handle(filename, loader, *args):
    def helper(fn):
        with fn(filename, 'rt') as f:
            return loader(f, *args)
    if filename == '-':
        return loader(sys.stdin, *args)
    _, ext = os.path.splitext(filename)
    if ext == '.bz2':
        return helper(bz2.open)
    elif ext == '.gz':
        return helper(gzip.open)
    else:
        return helper(open)

def load_snps(f, score_thresh):
    snps_raw = (line.split() for line in f)
    head = collections.defaultdict(list)
    for k, g in itertools.groupby(snps_raw, key=operator.itemgetter(0)):
        parsed = ((int(pos), name, float(score)) for _, pos, _, name, score in g)
        for p, n, s in parsed:
            if s >= score_thresh:
                head[k].append((p, n, s))
    return head

def load_regions(f):
    regions_raw = (line.split() for line in f)
    return {k: [(int(x[1]), int(x[2])) for x in g] for k, g in
            itertools.groupby(regions_raw, key=operator.itemgetter(0))}

def isect(snps, regions):
    region_iter = iter(regions)
    start, end = next(region_iter)
    for snp, name, score in snps:
        while end < snp:
            start, end = next(region_iter)
        if snp < start:
            continue
        if end > snp:
            yield name, score

def trial(snps, regions, offsets=None):
    if offsets is not None:
        regions = {k: [(s + o, e + o) for (s, e), o in zip(v, offsets)]
                   for k, v in regions.items()}
    overlaps = (x for k in snps for x in isect(snps[k], regions.get(k, [])))
    loci = itertools.groupby(sorted(overlaps, reverse=True), key=operator.itemgetter(0))
    return len([k for k, g in loci])

def moments(xs):
    running_mean = xs.pop(0)
    running_variance = 0
    for i, x in enumerate(xs):
        new_mean = running_mean + (x - running_mean) / (i + 2)
        running_variance += (x - running_mean) * (x - new_mean)
        running_mean = new_mean
    running_variance /= ntrials
    return running_mean, running_variance

def test(snps, regions, thresh):
    X = trial(snps, regions)
    offsets = ((-1 if random.random() > .5 else 1) * random.randrange(thresh)
               for _ in itertools.count())
    Y = [trial(snps, regions, offsets) for _ in range(ntrials)]
    mean, variance = moments(Y)
    _, p = scipy.stats.kstest(Y, scipy.stats.norm(mean, math.sqrt(variance)).cdf)
    exact_p = (1 + len([y for y in Y if y >= X])) / (1 + ntrials)
    return X, mean, variance, p, exact_p

if __name__ == '__main__':
    random.seed(0)
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
    if sys.argv[1] == sys.argv[2] == '-':
        print('error: SNPS and REGIONS cannot both be -')
        sys.exit(1)
    ntrials = int(sys.argv[3])
    if ntrials <= 0:
        print('error: N must be positive')
        sys.exit(1)
    thresh = int(sys.argv[4])
    if thresh <= 0:
        print('error: thresh must be positive')
        sys.exit(1)
    score_thresh = int(sys.argv[5])
    if score_thresh <= 0:
        print('error: score_thresh must be positive')
        sys.exit(1)
    snps = handle(sys.argv[1], load_snps, score_thresh)
    regions = handle(sys.argv[2], load_regions)
    result = test(snps, regions, thresh)
    print(' '.join('{:.3f}'.format(x) for x in result))
