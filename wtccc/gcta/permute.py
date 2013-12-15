import random
import sys

data = [line.split() for line in sys.stdin]
status = [d[2] for d in data]
random.shuffle(status)
for (fid, iid, _), s in zip(data, status):
    print(fid, iid, s)
    
