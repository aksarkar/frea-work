import sys

import scipy.stats

types = [str, int, int, int, int]
data = ([f(x) for f, x in zip(types, line.split())] for line in sys.stdin)
T = scipy.stats.fisher_exact
for name, fg_hits, fg_total, bg_hits, bg_total in data:
    print(name, *T([[fg_hits, fg_total - fg_hits],
                    [bg_hits, bg_total - bg_hits]],
                   alternative='greater'))
