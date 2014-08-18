"""Compute the number of LD partners for each SNP

Usage: python blocksize.py THRESH

Expects space-separated lines (snp1, snp2, r^2, d') on stdin

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import collections
import itertools
import sys

thresh = float(sys.argv[1])
data = (line.split() for line in sys.stdin)
high_pass = ((x, y) for x, y, r2, _ in data if float(r2) > thresh)
counts = collections.Counter(itertools.chain.from_iterable(high_pass))
for k in sorted(counts):
    print(k, counts[k])
