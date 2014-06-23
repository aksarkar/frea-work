"""Index 1KD LD files for quick lookups

Usage: python index.py LDFILE

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import operator
import sys

key = operator.itemgetter(*[int(x) for x in sys.argv[2:]])
with open(sys.argv[1], 'rb') as f:
    pos = 0
    line = f.readline()
    while line:
        print(str(key(line.split()), 'ascii'), pos)
        pos = f.tell()
        line = f.readline()
