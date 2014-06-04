import random
import sys

random.seed(sys.argv[1])
snps = [line.strip() for line in sys.stdin]
random.shuffle(snps)
for s in snps:
    print(s)
