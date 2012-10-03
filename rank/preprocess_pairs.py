import csv
import itertools
import math
import sys

import scipy.stats

r = csv.reader(sys.stdin)
xss = list(zip(*list(r)))
with open('pairs.csv', 'w') as f:
    print('disease1,disease2,rank1,rank2', file=f)
    for xs, ys in itertools.combinations(xss, 2):
        a, b = xs[0], ys[0]
        zs, ws = zip(*[(int(x), int(y)) for x, y in zip(xs[1:], ys[1:]) if
                       int(x) < 20000 or int(y) < 20000])
        for z, w in zip(zs, ws):
            print(a, b, z, w, sep=',', file=f)

with open('stats.csv', 'w') as f:
    print('disease1,disease2,cutoff,tt,tn,nt,nn,p,exp', file=f)
    cs = [10, 50, 100, 500, 1000, 2000, 3000, 4000, 5000, 10000, 15000, 20000]
    for c in cs:
        for xs, ys in itertools.combinations(xss, 2):
            a, b = xs[0], ys[0]
            zs, ws = zip(*[(int(x), int(y)) for x, y in zip(xs[1:], ys[1:]) if
                           int(x) < c or int(y) < c])
            tt = len([z for z, w in zip(zs, ws) if int(z) < c and int(w) < c])
            tn = len([z for z in zs if int(z) < c])
            nt = len([w for w in ws if int(w) < c])
            nn = len(xs) - 1 - nt
            exp = c * c / len(xs)
            _, p = scipy.stats.fisher_exact([[tt, tn], [nt, nn]])
            print(a, b, c, tt, tn, nt, nn, p, exp, sep=',', file=f)
