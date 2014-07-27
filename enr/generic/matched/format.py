import itertools
import operator
import sys

print('cluster', 'cells', 'z.score', 'empirical.p', sep='\t')
data = (line.split() for line in sys.stdin)
for k, g in itertools.groupby(data, key=operator.itemgetter(0,1,2)):
    cluster, z, p = k
    cells = sorted(set(x[3] for x in g))
    print(cluster, ','.join(cells), z, p, sep='\t')
