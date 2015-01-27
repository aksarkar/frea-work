"""Compute p-values for functional enrichment against resampled sets

Usage: python test.py TABLE THRESH NTRIALS

Expects BED file of SNPs with score = negative log p-value on stdin. TABLE is
gzipped space-separated (rsid, key1[, key2, ...]) tuples.

This implementation resamples SNPs exceeding THRESH and counts how many
overlaps occur over resampled sets.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import collections
import itertools
import gzip
import math
import random
import sys

import scipy.stats

def build_bins(snps, table, thresh):
    result = collections.defaultdict(list)
    counts = collections.defaultdict(int)
    for rsid, *key in table:
        k = tuple(key)
        result[k].append(rsid)
        if rsid in snps and snps[rsid][0] > thresh:
            counts[k] += 1
    return result, counts

I = itertools.chain.from_iterable
def match(bins, counts):
    return I((random.choice(bins[k]) for _ in range(v)) for k, v in
             counts.items() if k in bins)

def num_overlaps(rsids):
    return len([s for s in rsids if s in snps and snps[s][1]])

def moments(xs):
    running_mean = xs.pop(0)
    running_variance = 0
    for i, x in enumerate(xs):
        new_mean = running_mean + (x - running_mean) / (i + 2)
        running_variance += (x - running_mean) * (x - new_mean)
        running_mean = new_mean
    running_variance /= ntrials
    return running_mean, running_variance

def permutation_test(snps, bins, counts, thresh, ntrials):
    X = num_overlaps(s for s in snps if snps[s][0] > thresh)
    if X == 0:
        return 0, 0, 0, 1
    Y = [num_overlaps(match(bins, counts)) for _ in range(ntrials)]
    mean, variance = moments(Y)
    _, p = scipy.stats.kstest(Y, scipy.stats.norm(mean, math.sqrt(variance)).cdf)
    if p < .01:
        print('Warning: null distribution significantly different from normal (p={:.3f})'.format(p), file=sys.stderr)
    exact_p = (1 + len([y for y in Y if y >= X])) / (1 + ntrials)
    return X, mean, variance, exact_p

if __name__ == '__main__':
    random.seed(0)
    thresh = float(sys.argv[2])
    snps_raw = (line.split() for line in sys.stdin)
    snps = {row[3]: (float(row[4]), int(row[5])) for row in snps_raw}

    with gzip.open(sys.argv[1]) as f:
        data = (str(l, 'utf-8').split() for l in f)
        filter_ = (d for d in data if d[0] in snps)
        bins, counts = build_bins(snps, filter_, thresh)

    ntrials = int(sys.argv[3])
    result = permutation_test(snps, bins, counts, thresh, ntrials)
    print('{:.3f} {:.3f} {:.3f} {:.3f} {:.3g}'.format(thresh, *result))
