import random
import sys

random.seed(sys.argv[1])
data = [line.split() for line in sys.stdin]
cases = [(f, i) for f, i, c in data if c == '1']
controls = [(f, i) for f, i, c in data if c == '0']
for f, i in random.sample(cases, len(cases) // 2):
    print(f, i, 1)
for f, i in random.sample(controls, len(controls) // 2):
    print(f, i, 0)
