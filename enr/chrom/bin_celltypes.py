import csv
import itertools as it
import math
import sys

celltypes = ['GM12878','H1','HMEC','HSMM','HepG2','Huvec','K562','NHEK','NHLF']

r = csv.reader(sys.stdin)
head = next(r)
types = [str, float] + [int for _ in celltypes]
data = [[int(x) for x in xs[2:]] for xs in
        sorted(list(r), key=lambda row: float(row[1]))]

enhancers = 1
expected = [len([x for x in xs if x & enhancers]) / len(xs) for xs in
            zip(*data)]

binsize = int(sys.argv[1])
count = [0 for _ in celltypes]
w = csv.writer(sys.stdout)
w.writerow(['total', 'feature', 'count', 'expected'])
for i in range(0, len(data), binsize):
    count = [c + x for c, x in 
             zip(count, [len([x for x in xs if x & enhancers]) for
                         xs in zip(*data[i:i + binsize])])]
    for c, type_, e in zip(count, celltypes, expected):
        w.writerow([i + binsize, type_, c, int((i + binsize) * e)])
