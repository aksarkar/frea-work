import itertools
import operator
import sys

def mean(l):
    if not l:
        return 0
    else:
        return sum(l) / len(l)

binwidth = float(sys.argv[1])
types = [float, int, str]
data = sorted([[f(x) for f, x in zip(types, line.split())] for line in sys.stdin],
              key=operator.itemgetter(2, 0))
groups = itertools.groupby(data, key=operator.itemgetter(2))
for k, g in groups:
    curr = 0
    while curr < 1:
        curr += binwidth
        bin_ = [x for _, x, _ in itertools.takewhile(lambda x: x[0] < curr, g)]
        print(k, curr - binwidth / 2, mean(bin_))
