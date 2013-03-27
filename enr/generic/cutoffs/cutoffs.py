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
import itertools
import sys

import scipy.stats

markers = sys.argv[1]
feature = sys.argv[2]
celltype = sys.argv[3]
exclude = len(sys.argv) > 4

types = [float, int]
data = [tuple(g(x) for g, x in zip(types, line.split()))
        for line in sys.stdin]
if not data:
    raise ValueError('No data ({}, {}, {})'.format(markers, feature, celltype))
    
ns = list(itertools.takewhile(lambda m: m < len(data),
                              (100 << i for i in itertools.count())))
if exclude:
    ns.insert(0, 0)
    ns.append(2 * ns[-1])

num_overlaps = len([s for s, a in data if a])
F = scipy.stats.fisher_exact
if exclude:
    obs = [len([s for s, a in data[ra:min(rb, len(data))] if a])
           for ra, rb in zip(ns, ns[1:])]
    ms = [min(rb, len(data)) - ra for ra, rb in zip(ns, ns[1:])]
else:
    obs = [len([s for s, a in data[:n] if a]) for n in ns]
    ms = ns
tests = (F([[o, m - o], [num_overlaps - o, len(data) - m - num_overlaps + o]])
         for o, m in zip(obs, ms))

if exclude:
    ns.pop(0)
for n, (o, p) in zip(ns, tests):
    print(markers, feature, celltype, n, -math.log(p, 10), o)
