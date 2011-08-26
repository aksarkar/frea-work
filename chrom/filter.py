import csv
import sys

import bitvector

if __name__ == '__main__':
    cell_type = sys.argv[1]
    mask = bitvector.bitvector(int(x) for x in sys.argv[2][1:-1].split(','))
    ps, annots = zip(*sorted((float(row['p']), 
                              1 if int(row[cell_type]) & mask else 0)
                             for row in csv.DictReader(sys.stdin)))
    print(len(annots))
    for a in annots:
        print(a)
