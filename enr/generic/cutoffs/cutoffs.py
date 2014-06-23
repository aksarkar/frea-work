"""Compute enrichment at selected cutoffs

Usage: python cutoffs.py PHENOTYPE FEATURE CELLTYPE [EXCLUDE]

Reference: Maurano et al. "Systematic Localization of Common Disease-Associated
Variation in Regulatory DNA." Science. 2012. doi:10.1126/science.1222794

Expects whitespace-separated (score, binary annotation) sorted by decreasing
score on stdin. Writes space-separated (phenotype, feature, cell type, cutoff,
fold enrichment) on stdout.

If EXCLUDE is non-nil, exclude the SNPs meeting the previous cutoff (compute
enrichment for disjoint intervals of SNP ranks).

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import math
import collections
import itertools
import sys

import scipy.stats

phenotype = sys.argv[1]
feature = sys.argv[2]
celltype = sys.argv[3]
exclude = len(sys.argv) > 4

data = (line.split() for line in sys.stdin)
parsed = ((float(s), int(a)) for s, a in data)
overlaps = [0]
totals = []
cutoffs = []
breaks = (100 << i for i in itertools.count())
current_bin = next(breaks)
total_snps = 0
total_overlaps = 0
for i, (s, a) in enumerate(parsed):
    if i >= current_bin:
        current_bin = next(breaks)
        overlaps.append(0)
        totals.append(current_bin - i)
        cutoffs.append(s)
    if a:
        overlaps[-1] += 1
        total_overlaps += 1
    total_snps += 1
if not exclude:
    overlaps = itertools.accumulate(overlaps)
    totals = itertools.accumulate(totals)
for overlap_count, total_count, cutoff in zip(overlaps, totals, cutoffs):
    contingency = [[overlap_count, min(total_snps, total_count)],
                   [total_overlaps, total_snps]]
    odds_ratio, p = scipy.stats.fisher_exact(contingency, alternative='greater')
    logp = -math.log(p, 10) if p > 0 else 1000
    print('{} {} {} {:.3f} {:.3f} {:.3f}'.format(phenotype, feature, celltype,
                                                 cutoff, logp, odds_ratio))
