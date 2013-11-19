"""Generate PLINK map files matched to dosage files

Usage: python gmap.py RSIDS

Expects dosage on stdin. Writes map on stdout.

This is needed because FastLMM expects map entries to be in exactly the same
order as the dosage file. Hashing the map information takes less memory than
sorting the dosage file and then using join(1).

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import sys

with open(sys.argv[1]) as f:
    data = (line.split() for line in f)
    lookup = {row[1]: row for row in data}

for line in sys.stdin:
    print(' '.join(lookup[line.split()[0]]))
