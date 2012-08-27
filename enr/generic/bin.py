import csv
import sys

phenotype = sys.argv[1]
feature = sys.argv[2]
celltype = sys.argv[3]
binsize = int(sys.argv[4])
data = [int(x) for x in sys.stdin]
total = data.count(1)
w = csv.writer(sys.stdout)
w.writerow(['total', 'phenotype', 'celltype', 'feature', 'count', 'expected'])
count = 0
for i in range(0, len(data), binsize):
    count += data[i:i+binsize].count(1)
    w.writerow([i + binsize, phenotype, celltype, feature, count,
                int((i + binsize) * total / len(data))])
