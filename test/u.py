import sys

import scipy.stats

types = [float, int]
data = [tuple(g(x) for g, x in zip(types, line.split())) for line in sys.stdin]
if len(set(c for p, c in data)) == 1:
    print(0)
    print('Warning: all markers in one class', file=sys.stderr)
else:
    u, prob = scipy.stats.mannwhitneyu([p for p, c in data if c == 1],
                                       [p for p, c in data if c == 0])
    print('{:6e} {:6e}'.format(u, prob))
