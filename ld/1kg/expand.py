"""Expand SNPs according to 1KG LD

Usage: python expand.py TAGS

Expects id1, id2 pairs on stdin. Writes pairs on stdout

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""

import sys

with open(sys.argv[1]) as f:
    tags = set(line.strip() for line in f)

data = (line.split() for line in sys.stdin)
for row in data:
    if row[0] in tags:
        print(row[0], row[1])
    elif row[1] in tags:
        print(row[1], row[0])
