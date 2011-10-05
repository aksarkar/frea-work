import re
import sys

import pymysql

ld_query = 'select ld from ceu_ld where id={}'
pos_query = 'select hg19chr, hg19pos from hg19liftover where id={}'
patt = re.compile('[,;]')

def parse(s, rsid):
    t = patt.split(s)
    if len(t) == 1:
        t = []
    t.extend([rsid, 1])
    # print(t, file=sys.stderr)
    return zip(t[::2], t[1::2])

def lookup_pos(rsid, cur):
    cur.execute(pos_query.format(conn.escape(rsid)))
    if not cur.rowcount:
        return []
    else:
        return [rsid] + [int(x) for x in cur.fetchone()]

if __name__ == '__main__':
    conn = pymysql.connect(db='HaploReg', host='calcium', user='HaploReg',
                           passwd='HaploReg')
    cur = conn.cursor()
    thresh = float(sys.argv[1])
    for line in sys.stdin:
        rsid, pv = line.split()
        cur.execute(ld_query.format(conn.escape(rsid)))
        ld = cur.fetchone()[0] if cur.rowcount else ''
        ps = [lookup_pos(marker, cur) for marker, ld in parse(ld, rsid)
              if float(ld) >= thresh]
        for p in ps:
            if p:
                print('chr{}'.format(p[1]), p[2], p[2] + 1, rsid, pv, sep='\t')
    conn.close()
