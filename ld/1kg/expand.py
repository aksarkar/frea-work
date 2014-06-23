import csv
import sys

with open(sys.argv[1]) as f:
    tags = set(line.strip() for line in f)

for row in csv.reader(sys.stdin, delimiter=' '):
    if row[0] in tags:
        print(row[0], row[1])
    elif row[1] in tags:
        print(row[1], row[0])
