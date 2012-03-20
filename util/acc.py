#!/usr/bin/env python3
import collections
import csv
import sys

def entries(f):
    return {(int(row[0]), row[1]): int(row[2]) for row in csv.reader(open(f))}

c = collections.Counter(entries(sys.argv[1]))
c.update(entries(sys.argv[2]))
for (bin, annot), count in sorted(c.items()):
    print(bin, annot, count, sep=',')
