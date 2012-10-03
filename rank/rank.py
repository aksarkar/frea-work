import csv
import sys

xs = [float(row[-1]) for row in csv.reader(sys.stdin, delimiter='\t')]
ranks = [y[0] for y in
         sorted(enumerate(sorted((x, i) for i, x in enumerate(xs))),
                key=lambda x: x[1][1])]
for r in ranks:
    print(r)
