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
import itertools
import sys

markers = sys.argv[1]
feature = sys.argv[2]
celltype = sys.argv[3]
exclude = len(sys.argv) > 4

types = [float, int]
data = [tuple(g(x) for g, x in zip(types, line.split()))
        for line in sys.stdin]
if not data:
    raise ValueError('No data ({}, {}, {})'.format(markers, feature, celltype))
    
ns = list(itertools.takewhile(lambda m: m < 2 * len(data),
                              (100 << i for i in itertools.count())))
if exclude:
    ns.insert(0, 1)
cutoffs = [data[n - 1][0] if n - 1 < len(data) else data[-1][0] for n in ns]

positives = [s for s, a in data if a]
if exclude:
    obs = [len([s for s in positives if sb <= s <= sa]) / (rb - ra)
           for sa, sb, ra, rb in zip(cutoffs, cutoffs[1:], ns, ns[1:])]
else:
    obs = [len([s for s in positives if c <= s]) / min(n, len(data))
           for c, n in zip(cutoffs, ns)]

exp = (1 + len(positives)) / (1 + len(data))
folds = (o / exp for o in obs)
if exclude:
    ns.pop(0)
    cutoffs.pop(0)
for n, c, f in zip(ns, cutoffs, folds):
    print(markers, feature, celltype, n, c, f)
