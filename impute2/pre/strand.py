import sys

ambiguous = ('AT', 'TA', 'CG', 'GC')

data = (line.split() for line in sys.stdin)
for row in data:
    rsid, _, _, a0, a1, *_ = row
    if a0 + a1 in ambiguous:
        print(rsid)
