import csv
import math
import sys

if __name__ == '__main__':
    filter_comment = (x for x in sys.stdin if x.strip() and x[0] != '#')
    reader = csv.DictReader(filter_comment, delimiter='\t')
    filter_missing = (row for row in reader if row['P-value'] and
                      row['Chr Position'] and row['Chr ID'])
    for row in filter_missing:
        print('chr{}'.format(row['Chr ID']),
              row['Chr Position'],
              int(row['Chr Position']) + 1,
              row['SNP ID'],
              -math.log10(float(row['P-value'])), sep='\t')
