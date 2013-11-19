"""Join dosage files across cohorts

Usage: python join.py MAP DOSAGE

MAP is the corresponding map file (SNPs not required to be in the same
order). Expects one dosage file on stdin (useful for constructing pipelines to
join multiple cohorts).

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import csv
import sys

with open(sys.argv[1]) as f:
    pos = {row[1]: int(row[3]) for row in csv.reader(f, delimiter=' ')}

with open(sys.argv[2]) as f:
    r1 = (l.split() for l in sys.stdin)
    r2 = (line.split() for line in f)
    a = next(r1)
    b = next(r2)
    while True:
        try:
            while b[0] not in pos:
                b = next(r2)
            while pos.get(a[0], 0) < pos[b[0]]:
                a = next(r1)
            while pos[a[0]] > pos.get(b[0], 0):
                b = next(r2)
            if pos[a[0]] == pos[b[0]]:
                print(*(a + b[3:]))
            a = next(r1)
        except StopIteration:
            sys.exit(0)
