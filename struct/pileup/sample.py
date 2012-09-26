import csv
import random
import sys

n = int(sys.argv[1])
for x in random.sample(list(csv.reader(sys.stdin, delimiter='\t')), n):
    x[3] = x[1]
    print(*x, sep='\t')
