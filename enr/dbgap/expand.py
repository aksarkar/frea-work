import collections
import csv
import sys

thresh = float(sys.argv[2])
tags = {x[0]: float(x[1]) for x in csv.reader(sys.stdin, delimiter=' ')}

snps = {t: set(t) for t in tags}
with open(sys.argv[1]) as f:
    for row in csv.reader(f, delimiter=' '):
        if row[0] in snps:
            snps[row[0]].add(row[1])
        elif row[1] in snps:
            snps[row[1]].add(row[0])
for s in snps:
    for t in snps[s]:
        print(s, t, tags[s])
