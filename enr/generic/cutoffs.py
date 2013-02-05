"""Compute enrichment at selected cutoffs

Usage: python cutoffs.py PHENOTYPE FEATURE CELLTYPE [EXCLUDE]

Reference: Maurano et al. "Systematic Localization of Common Disease-Associated
Variation in Regulatory DNA." Science. 2012. doi:10.1126/science.1222794

Expects whitespace-separated (p-value, binary annotation) on stdin. Writes
space-separated (phenotype, feature, cell type, cutoff, fold enrichment) on
stdout.

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
data = sorted(tuple(g(x) for g, x in zip(types, line.split()))
              for line in sys.stdin)

fracs = itertools.takewhile(lambda x: 2 ** x < len(data), itertools.count())
ms = (len(data) // (2 ** i) for i in fracs)
ns = list(itertools.takewhile(lambda n: n > 100, ms))
if exclude:
    ns.append(1)
cutoffs = [data[n - 1][0] for n in ns]

positives = [p for p, x in data if x]
if exclude:
    obs = ((len([p for p in positives if lower_p <= p < upper_p]) /
            (upper_rank - lower_rank))
           for lower_p, upper_p, lower_rank, upper_rank in
           zip(cutoffs[1:], cutoffs, ns[1:], ns))
else:
    obs = (len([p for p in positives if p < c]) / n for c, n in
           zip(cutoffs, ns))

exp = len(positives) / len(data)
folds = (o / exp for o in obs)
for n, c, f in zip(ns, cutoffs, folds):
    print(markers, feature, celltype, n, c, f)
