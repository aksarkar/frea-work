"""Compute hypergeometric p-values for functional enrichment with MAF matched
negative sets

Usage: python fisher.py TABLE THRESH NTRIALS

Expects BED file of SNPs with score = negative log p-value on stdin. TABLE is
gzipped space-separated (rsid, key1[, key2, ...]) tuples.

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

I = itertools.chain.from_iterable

def build_bins(snps, table):
    overlap_bins = collections.defaultdict(list)
    nonoverlap_bins = collections.defaultdict(list)
    for rsid, *key in table:
        k = tuple(key)
        if snps[rsid][1]:
            overlap_bins[k].append(rsid)
        else:
            nonoverlap_bins[k].append(rsid)
    return overlap_bins, nonoverlap_bins

def match(overlap_counts, nonoverlap_bins):
    return I((random.choice(nonoverlap_bins[k]) for _ in range(v)) for k, v in
             overlap_counts.items() if k in nonoverlap_bins)

def num_top_snps(rsids, thresh):
    result = 0
    for s in rsids:
        if snps[s][0] > thresh:
            result += 1
    return result

def permutation_test(overlap_bins, nonoverlap_bins, thresh, ntrials):
    X = num_top_snps(I(overlap_bins.values()), thresh)
    success_count = 0
    running_mean = X
    running_variance = 0
    overlap_counts = {k: len(overlap_bins[k]) for k in overlap_bins}
    for i in range(ntrials):
        Y = num_top_snps(match(overlap_counts, nonoverlap_bins), thresh)
        new_mean = running_mean + (Y - running_mean) / (i + 2)
        running_variance += (Y - running_mean) * (Y - new_mean)
        running_mean = new_mean
        if Y >= X:
            success_count += 1
    z = (X - running_mean) / math.sqrt(running_variance / ntrials)
    fold = X / running_mean
    asymptotic_p = 1 - scipy.stats.norm.cdf(z)
    exact_p = success_count / ntrials
    return z, fold, asymptotic_p, exact_p

if __name__ == '__main__':
    random.seed(0)
    snps_raw = (line.split() for line in sys.stdin)
    snps = {row[3]: (float(row[4]), int(row[5])) for row in snps_raw}

    with gzip.open(sys.argv[1]) as f:
        data = (str(l, 'utf-8').split() for l in f)
        filter_ = (d for d in data if d[0] in snps)
        overlap_bins, nonoverlap_bins = build_bins(snps, filter_)

    thresh = float(sys.argv[2])
    ntrials = int(sys.argv[3])
    result = permutation_test(overlap_bins, nonoverlap_bins, thresh, ntrials)
    print(thresh, *result)
