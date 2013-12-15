import csv
import sys

import scipy.stats

with open(sys.argv[1]) as f:
    enh_counts = {row[0]: int(row[1]) for row in csv.reader(f, delimiter=' ')}

with open(sys.argv[2]) as f:
    motif_counts = {row[0]: int(row[1]) for row in csv.reader(f, delimiter=' ')}

with open(sys.argv[3]) as f:
    isect_counts = {tuple(row[:2]): int(row[2]) for row in csv.reader(f, delimiter=' ')}

with open(sys.argv[4]) as f:
    domain_counts = {row[0]: int(row[1]) for row in csv.reader(f, delimiter=' ')}

for (cell, motif), n in isect_counts.items():
    m = domain_counts[motif] - enh_counts[cell] - motif_counts[motif] + n
    odds, p = scipy.stats.fisher_exact([[n, enh_counts[cell] - n],
                                        [motif_counts[motif] - n, m]])
    print(cell, motif, odds, p)
