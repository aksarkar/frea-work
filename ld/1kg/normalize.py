import operator
import sys

for line in sys.stdin:
    id_, inld = line.split()
    snps = (i.split(',') for i in inld.split(';'))
    types = (str, float, float)
    casted = ([f(x) for f, x in zip(types, s)] for s in snps)
    for i in sorted(casted, key=operator.itemgetter(2), reverse=True):
        print(id_, *i)
