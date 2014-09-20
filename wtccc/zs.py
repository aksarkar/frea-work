import csv
import math
import sys

data = csv.DictReader(sys.stdin, delimiter='\t')
for row in data:
    p = float(row['frequentist_add'])
    if p > 0 and p < 1:
        print(row['id'], -math.log10(p))
