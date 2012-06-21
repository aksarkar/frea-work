"""Logistic regression for WTCCC T1D

Computes ROC curves for L1-penalized logistic regression classifiers with
stratified 20-fold cross validation.

Usage: python2 classifier.py GENOTYPES

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
from __future__ import print_function
from __future__ import division

import itertools
import sys

import numpy
import sklearn.cross_validation as cv
import sklearn.linear_model as lm
import sklearn.metrics as m

A = numpy.array
with open(sys.argv[1]) as f:
    rows = [[int(x) for x in xs.split()[1:]] for xs in f]
    X = A(zip(*rows))
y = A([1 for _ in xrange(2000)] + [0 for _ in xrange(3004)])
model = lm.LogisticRegression(penalty='l1')
for i, (train, test) in enumerate(cv.StratifiedKFold(y, 20)):
    model.fit(X[train], y[train])
    ps = model.predict_proba(X[test])[:,1]
    with open('roc.{:02d}'.format(i), 'w') as f:
        for t in sorted(set(ps)):
            y_ = [1 if p > t else 0 for p in ps]
            tp, fn, fp, tn = m.confusion_matrix(y[test], y_).flatten()
            tpr = tp / (tp + fn)
            fpr = fp / (fp + tn)
            print(fpr, tpr, m.fbeta_score(y[test], y_, .5), file=f)
