import csv
import operator
import sys

with open(sys.argv[1]) as f:
    lookup = {row[0]: row[1] for row in csv.reader(f, delimiter='\t')}

types = [int, str, str, str, float]
data = [[f(x) for f, x in zip(types, r)] for r in csv.reader(sys.stdin)]
point = [r for r in data if r[0] == 20000]
pass_ = set([r[1] for r in sorted(point, reverse=True, key=operator.itemgetter(-1))][:100])
w = csv.writer(sys.stdout, quoting=csv.QUOTE_MINIMAL)
for r in data:
    if r[1] in pass_:
        r[1] = lookup[r[1]]
        w.writerow(r)
