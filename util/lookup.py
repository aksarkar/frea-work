import sys

import pymysql

query = {'hg18': 'select hg18chr, hg18pos from snp where id="{}";',
         'hg19': 'select hg19chr, hg19pos from hg19liftover where id="{}";'
         }[sys.argv[1]]

with pymysql.connect(db='HaploReg', host='calcium', user='HaploReg', passwd='HaploReg') as cur:
    for line in sys.stdin:
        cur.execute(query.format(line.strip()))
        for row in cur:
            print('chr{}'.format(row[0]), int(row[1]) - 1, int(row[1]), line.strip(), sep='\t')
