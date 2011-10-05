import csv
import sys

with open(sys.argv[1]) as f:
    top = set(row[0] for row in csv.reader(f, delimiter='\t'))
with open(sys.argv[2]) as f:
    w = csv.writer(sys.stdout, delimiter='\t')
    w.writerows(row for row in csv.reader(f, delimiter='\t')
                if row[3] not in top)
