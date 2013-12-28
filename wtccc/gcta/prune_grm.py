import collections
import sys

with open(sys.argv[1]) as f:
    samples = [line.strip() for line in f]

pairs = collections.defaultdict(set)
data = (line.split() for line in sys.stdin)
for i, j in data:
    pairs[i].add(j)
    pairs[j].add(i)
while pairs:
    best = sorted(pairs.keys(), key=lambda x: len(pairs[x]), reverse=True)[0]
    print(samples[int(best) - 1])
    pairs.pop(best)
    for ind in pairs:
        if best in pairs[ind]:
            pairs[ind].remove(best)
    pairs = {k: v for k, v in pairs.items() if v}
