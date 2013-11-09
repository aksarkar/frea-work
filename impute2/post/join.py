import csv
import sys

with open(sys.argv[1]) as f:
    pos = {row[1]: int(row[3]) for row in csv.reader(f, delimiter=' ')}

with open(sys.argv[2]) as f:
    r1 = csv.reader(sys.stdin, delimiter=' ')
    r2 = csv.reader(f, delimiter=' ')
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
