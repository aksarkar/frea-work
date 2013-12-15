import random
import sys

data = [line.split() for line in sys.stdin]
cases = [(f, i) for f, i, c in data if c == '1']
controls = [(f, i) for f, i, c in data if c == '0']
k = int(sys.argv[1])
for i in range(k):
    with open('{}.txt'.format(i), 'w') as out:
        for f, i in random.sample(cases, len(cases) // 2):
            print(f, i, 1, file=out)
        for f, i in random.sample(controls, len(controls) // 2):
            print(f, i, 0, file=out)
