import csv
import os
import sys

with open('{}/annot'.format(os.getenv('PT_WORK'))) as f:
    lookup = set(row['rsid'] for row in csv.DictReader(f) if 
                 row['cell_type'] == sys.argv[1])
ps, annots = zip(*sorted((float(row['p']), 
                          1 if row['rsid'] in lookup else 0)
                         for row in csv.DictReader(sys.stdin)))
print(len(annots))
for a in annots:
    print(a)
