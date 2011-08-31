import csv
import gzip
import sys

lookup = set(row['rsid'] for row in csv.DictReader(sys.stdin) if 
             row['cell_type'] == sys.argv[1])
with gzip.open(sys.argv[2]) as f:
    decode = (str(x, 'ascii') for x in f)
    filter_comment = (x for x in decode if x.strip() and x[0] != '#')
    r = csv.DictReader(filter_comment, delimiter='\t')
    ps, annots = zip(*sorted((float(row['P-value']), 
                              1 if row['SNP ID'] in lookup else 0)
                             for row in r))
print(len(annots))
for a in annots:
    print(a)
