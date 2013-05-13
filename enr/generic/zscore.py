import csv
import itertools
import operator
import sys

# import scipy.stats

print(sys.stdin.encoding)

# key = lambda x: x[:4]
# data = (line.split(',') for line in sys.stdin)
# for k, g in itertools.groupby(sorted(data, key=key), key=key):
#     xs = list(g)
#     fs = [float(x[4]) / float(x[5]) for x in xs]
#     zs = scipy.stats.zscore(fs)
#     for x, z in zip(xs, zs):
#         if len(x) == 6:
#             print(*k, sep=',', end=',')
#             print(z)
