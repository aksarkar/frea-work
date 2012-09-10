import itertools
import math
import sys

import scipy.stats

markers = sys.argv[1]
feature = sys.argv[2]
types = [float, int]
data = sorted(tuple(g(x) for g, x in zip(types, line.split()))
              for line in sys.stdin)
positives = [p for p, x in data if x]
denom = len(positives) / len(data)
fracs = itertools.takewhile(lambda x: 2 ** x < len(data), itertools.count(1))
ms = (len(data) // (2 ** i) for i in fracs)
ns = list(itertools.takewhile(lambda n: n > 100, ms))
cutoffs = [data[n][0] for n in ns]
folds = (len([p for p in positives if p < c]) / n / denom
         for c, n in zip(cutoffs, ns))
for c, f in zip(cutoffs, folds):
    print(markers, feature, c, f)
