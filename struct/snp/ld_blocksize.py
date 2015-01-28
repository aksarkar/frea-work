"""Count LD partners

Usage: python ld_blocksize.py

Expects (id1, id2, r, d') on stdin. Writes counts on stdout

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import collections
import itertools
import sys

data = (line.split()[:2] for line in sys.stdin)
ids = itertools.chain.from_iterable(data)
counts = collections.Counter(ids)
for k in sorted(counts):
    print(k, counts[k])
