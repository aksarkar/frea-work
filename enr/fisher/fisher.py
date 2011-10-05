#!/bin/env python
import sys

import scipy.stats

nx = int(next(sys.stdin))
ny = int(next(sys.stdin))

trans = {'Strong_Enhancer': 'strong',
         'Weak_Enhancer': 'weak',
         'Enhancer': 'all'}

print('cell_type,state,p,fold')
for line in sys.stdin:
    ctype, enh = line.strip().split()
    xa = int(next(sys.stdin))
    xb = nx - xa
    ya = int(next(sys.stdin))
    yb = ny - ya
    exp = nx * (xa + ya) / (nx + ny)
    print(ctype, trans[enh],
          scipy.stats.fisher_exact([[xa, xb], [ya, yb]],
                                   alternative='greater')[1],
          xa / exp, sep=',')
