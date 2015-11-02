import itertools
import operator
import random
import sys

random.seed(sys.argv[1])
data = sorted([line.split() for line in sys.stdin], key=operator.itemgetter(2))
cohorts = {k: [x[:2] for x in g] for k, g in
           itertools.groupby(data, key=operator.itemgetter(2))}
for k in cohorts:
    random.shuffle(cohorts[k])
for k in cohorts:
    for i, x in enumerate(cohorts[k]):
        if i % 2:
            print(*x)
