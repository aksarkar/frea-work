"""Filter IMPUTE2 format dosages based on rsid

Usage: python filter_dosages.py

Expects (rsid, file, offset) pairs on stdin

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import itertools
import operator
import sys

index = (line.split() for line in sys.stdin)
for k, g in itertools.groupby(index, key=operator.itemgetter(1)):
    with open(k, 'rb') as f:
        for _, _, o in g:
            f.seek(int(o))
            print(str(f.readline(), 'utf-8'), end='')
