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
    snps = set(line.strip() for line in f)
    print('loaded {} SNPs'.format(len(snps)), file=sys.stderr)
with open(sys.argv[2]) as f:
    ref_ld = (line.split() for line in f)
    ld = collections.defaultdict(set)
    for k, v, r, d in ref_ld:
        if float(r) > thresh and k in snps and v in snps:
            ld[k].add(v)
            ld[k].add(k)
            ld[v].add(k)
            ld[v].add(v)
    print('found {} pairs in LD'.format(len(ld)), file=sys.stderr)

num_tagged = lambda x: len(ld[x])
n = 0
while ld:
    tag = sorted(ld.keys(), key=num_tagged, reverse=True)[0]
    n += 1
    assert n < len(snps)
    prune = ld[tag]
    for snp in prune:
        print(tag, snp)
        ld.pop(snp)
    for snp in ld:
        ld[snp] = ld[snp].difference(prune)
