import csv
import itertools
import sys

window_size = 5000000
data = (line.split() for line in sys.stdin)
next(data)
pos = (int(row[2]) for row in data)
windows = (k for k, g in
           itertools.groupby(pos, key=lambda x: x // window_size)
           if list(g))
for i, k in enumerate(windows):
    print('{:02d} {} {}'.format(i, k * window_size, (k + 1) * window_size))
