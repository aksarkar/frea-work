"""Prune related pairs of individuals from the GRM

Usage: python prune_grm.py FAM

Based on the implementation in GCTA 1.20. Expects (ID1, ID2, N,
relatedness) tuples on stdin. Writes individuals to be excluded on
stdout.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import itertools
import sys

with open(sys.argv[1]) as f:
    samples = [line.strip() for line in f]

data = (line.split() for line in sys.stdin)
grm = ((int(a), int(b), float(k)) for a, b, _, k in data)
A, B = [list(x) for x in zip(*[(a, b) for a, b, k in grm
                               if a > b and k > .025])]
C = {k: len(list(g)) for k, g in itertools.groupby(sorted(A + B))}
for i in range(len(A)):
    if C[A[i]] < C[B[i]]:
        A[i], B[i] = B[i], A[i]
prune = (k for k, _ in itertools.groupby(sorted(A)))
for p in prune:
    print(samples[p - 1])
