import collections
import sys

data = (line.split() for line in sys.stdin)
ld = collections.defaultdict(set)
for a, b in data:
    ld[a].add(b)
    ld[b].add(a)

with open(sys.argv[1]) as f:
    ranking = [line.strip() for line in f]

tags = set()
tagged = set()
for i, x in enumerate(ranking):
    if x not in tagged:
        tags.add(x)
        tagged |= ld[x]
    if i > 0 and not i % 1000:
        print(i, len(tags))
