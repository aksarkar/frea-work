"""Prune SNPs according to LD

Usage: python prune.py SNPS LD [THRESH]

Iteratively pick the SNP which tags the most other SNPs to add to our list,
then remove it and its LD partners until no pairs of SNPs in LD remains.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import collections
import sys

thresh = float(sys.argv[3]) if len(sys.argv) > 3 else 0
with open(sys.argv[1]) as f:
    types = [str, int, int, str, float]
    scores = (line.split() for line in f)
    parsed = ([f(x) for f, x in zip(types, snp)] for snp in scores)
    snps = {name: score for _, _, _, name, score in parsed}

with open(sys.argv[2]) as f:
    ref_ld = (line.split() for line in f)
    ld = collections.defaultdict(list)
    for k, v, r, d in ref_ld:
        if float(r) > thresh and k in snps and v in snps:
            ld[k].append((v, r))
            ld[v].append((k, r))

num_tagged = lambda x: len(ld[x])
n = 0
while ld:
    tag = sorted(ld.keys(), key=num_tagged, reverse=True)[0]
    n += 1
    assert n <= len(snps)
    prune = ld[tag]
    for snp, r in prune:
        if snp in ld:
            print(tag, snp, r)
            ld.pop(snp)
