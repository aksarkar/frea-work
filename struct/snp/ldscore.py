"""Compute LD score

Usage: python ldscore.py

Expects (id1, id2, r, d') on stdin. Writes (id, score) pairs on stdout

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import collections
import itertools
import sys

data = (line.split() for line in sys.stdin)
scores = collections.defaultdict(int)
for a, b, r, _ in data:
    r = float(r)
    scores[a] += r
    scores[b] += r
for k in sorted(scores):
    print(k, scores[k])
