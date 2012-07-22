"""Index 1KD LD files for quick lookups

Usage: python index.py LDFILE

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import sys

with open(sys.argv[1], 'rb') as f:
    pos = 0
    line = f.readline()
    while line:
        pos = f.tell() - len(line)
        print(str(line.split()[0], 'ascii'), pos)
        line = f.readline()
