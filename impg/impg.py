import collections
import csv
import sys

with open(sys.argv[1]) as f:
    r = csv.DictReader(f, delimiter='\t')
    typed_zscores = {row['SNP']: float(row['SNPWeight']) / float(row['SNPWeightSE'])
                     for row in r}

for beta_file in sys.argv[2:]:
    with open(beta_file) as f:
        data = (line.split() for line in f)
        header = next(data)[2:-1]
        for rsid, _, *rest in data:
            if rsid in typed_zscores:
                imputed_zscore = typed_zscores[rsid]
                r2pred = 1
            else:
                parsed = (float(x) for x in rest)
                imputed_zscore = sum(typed_zscores[typed_rsid] * b
                                     for typed_rsid, b in
                                     zip(header, parsed))
                r2pred = next(parsed)
            print('{} {:.3f} {:.3f}'.format(rsid, imputed_zscore, r2pred))
            
