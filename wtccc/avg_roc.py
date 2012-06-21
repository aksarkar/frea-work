"""Compute a threshold average ROC curve

Usage: python3 avg_roc.py NSAMPLES LABEL ROC1 [...]

Expects ROC curves as CSV with entries (fpr, tpr, score) sorted by fpr. Writes
out CSV with entries (fpr, tpr, label).

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import itertools
import functools
import math
import operator
import sys

def load(filename):
    with open(filename) as f:
        return sorted(([float(x) for x in line.split()] for line in f),
                      key=operator.itemgetter(2))

def mean(xs):
    return sum(xs) / len(xs)

def centroid(xs, ys, zs):
    return mean(xs), mean(ys), mean(zs)

def at(rocs, thresh):
    def helper(r, t):
        return max((s, i) for i, (_, _, s) in enumerate(r) if s <= t)[1]
    indexes = (helper(r, thresh) for r in rocs)
    return [r[i] for r, i in zip(rocs, indexes)]

def area(p, q, r):
    px, py, *_ = p
    qx, qy, *_ = q
    rx, ry, *_ = r
    return 5 * (-qx * py + rx * py + px * qy - rx * qy - px * ry + qx * ry)

def hull(ps):
    h, ps = ps[:2], ps[2:]
    while ps:
        if len(h) > 1:
            if area(h[-2], h[-1], ps[0]) < 0:
                h.append(ps.pop(0))
            else:
                h.pop()
        else:
            h.append(ps.pop(0))
    return h

nsamples = int(sys.argv[1])
label = sys.argv[2]
rocs = [load(a) for a in sys.argv[3:]]
scores = sorted(s for x, y, s in itertools.chain(*rocs))
thresholds = (scores[i] for i in range(0, len(scores), len(scores) // nsamples))
samples = sorted(centroid(*list(zip(*at(rocs, t)))) for t in thresholds)
for x, y, z in hull(samples):
    print(x, y, label, sep=',')
