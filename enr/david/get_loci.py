from __future__ import print_function

import itertools
import operator
import signal
import sys

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

thresh = float(sys.argv[1])

data = (line.split() for line in sys.stdin)
types = [str, int, int, str, float]
parsed = ([f(x) for f, x in zip(types, row)] for row in data)
for k, g in itertools.groupby(parsed, key=operator.itemgetter(3)):
    snps = list(g)
    tag = max(x[-1] for x in snps) 
    if thresh < tag < 8:
        for x in snps:
            print(*x, sep='\t')
