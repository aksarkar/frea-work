from __future__ import print_function

import sklearn.cross_validation as cv

y = [1 for _ in range(2000)] + [0 for _ in range(3004)]
c = cv.StratifiedKFold(y, k=20)
with open('train', 'w') as f, open('test', 'w') as g:
    for i, (train, test) in enumerate(c):
        for t in train:
            print(i, t, sep=',', file=f)
        for t in test:
            print(i, t, sep=',', file=g)
