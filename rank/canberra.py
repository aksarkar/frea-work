import sys

import scipy.stats
import scipy.spatial.distance

data = [line.split() for line in sys.stdin]
parsed = ((float(row[3]), float(row[4])) for row in data)
print(scipy.stats.spearmanr(*zip(*parsed)))
# print(scipy.spatial.distance.canberra(*zip(*parsed)) / len(data))
