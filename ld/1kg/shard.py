import gzip
import itertools
import operator
import sys

with gzip.open(sys.argv[1], 'rt') as f:
    data = (line.split() for line in f)
    # for k, g in itertools.groupby(data, key=lambda x: x[0].split('|')[1]):
    for k, g in itertools.groupby(data, key=operator.itemgetter(0)):
        with gzip.open('{}.txt.gz'.format(k), 'wt') as out:
            for row in g:
                print(*row, file=out)
    
