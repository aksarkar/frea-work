from __future__ import print_function

import sys

import scipy.stats

types = [float, int]
data = [tuple(g(x) for g, x in zip(types, line.split())) for line in sys.stdin]
d, prob = scipy.stats.ks_2samp([p for p, c in data if c == 1],
                               [p for p, c in data if c == 0])
print('{:6e}'.format(prob))
