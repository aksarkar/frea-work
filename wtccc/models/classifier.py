"""Logistic regression for WTCCC T1D

Performs L1-penalized logistic regression with 20-fold stratified
cross-validation. Outputs model coefficients, probabilities of positive class,
individual ROC curves, and threshold averaged ROC curve.

Usage: python2 classifier.py GENOTYPES TRAIN TEST

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
from __future__ import print_function
from __future__ import division

import csv
import functools
import itertools
import operator
import sys

import numpy
import sklearn.cross_validation as cv
import sklearn.linear_model as lm
import sklearn.metrics as m

def roc(ps, y, weight=.5):
    """Yield points (fpr, tpr, beta) on the ROC curve"""
    yield 0, 0, 0
    for t in sorted(set(ps)):
        y_ = [1 if p > t else 0 for p in ps]
        tp, fn, fp, tn = m.confusion_matrix(y, y_).flatten()
        tpr = tp / (tp + fn)
        fpr = fp / (fp + tn)
        yield fpr, tpr, m.fbeta_score(y, y_, weight)

def mean(xs):
    """Return the mean of list xs"""
    return sum(xs) / len(xs)

def centroid(xs, ys, zs):
    """Return the centroid of points (x, y, z)"""
    return mean(xs), mean(ys), mean(zs)

def at(rocs, thresh):
    """Return the greatest point with score less than thresh for each ROC
    curve"""
    def helper(r, t):
        return max((s, i) for i, (_, _, s) in enumerate(r) if s <= t)[1]
    indexes = (helper(r, thresh) for r in rocs)
    return [r[i] for r, i in zip(rocs, indexes)]

def avg_roc(rocs, nsamples=200):
    """Yield points (fpr, tpr, beta) on the treshold averaged ROC curve"""
    scores = sorted(s for x, y, s in itertools.chain(*rocs))
    thresholds = (scores[i] for i in range(0, len(scores),
                                           len(scores) // nsamples))
    yield 0, 0, 0
    for t in thresholds:
        yield centroid(*list(zip(*at(rocs, t))))
    yield 1, 1, 0

def area(p, q, r):
    """Return the signed area of the triangle (p, q, r)"""
    px, py = p[:2]
    qx, qy = q[:2]
    rx, ry = r[:2]
    return 5 * (-qx * py + rx * py + px * qy - rx * qy - px * ry + qx * ry)

def hull(ps):
    """Return the convex hull of list of points ps.

    Assume points are sorted by x, monotonically increasing in y.

    """
    h, ps = ps[:2], ps[1:]
    while ps:
        if len(h) > 1:
            if area(h[-2], h[-1], ps[0]) < 0:
                h.append(ps.pop(0))
            else:
                h.pop()
        else:
            h.append(ps.pop(0))
    return h

def run(model, X, y, train, test):
    """Run one iteration of cross-validation.

    Return model coefficients, predictions, and ROC curve

    """
    model.fit(X[train], y[train])
    ps = model.predict_proba(X[test])[:,1]
    return model.coef_, ps

def load_set(f):
    return [[int(x[1]) for x in xs] for k, xs in
            itertools.groupby(csv.reader(f), key=operator.itemgetter(0))]

if __name__ == '__main__':
    A = numpy.array
    with open(sys.argv[1]) as f:
        X = A(zip(*[[int(x) for x in xs.split()[1:]] for xs in f]))
    y = A([1 for _ in xrange(2000)] + [0 for _ in xrange(3004)])
    model = lm.LogisticRegression(penalty='l1')
    run_ = functools.partial(run, model, X, y)
    with open(sys.argv[2]) as f:
        train_sets = load_set(f)
    with open(sys.argv[3]) as f:
        test_sets = load_set(f)
    css, pss = zip(*[run_(train, test) for train, test in
                     zip(train_sets, test_sets)])
    with open('models.csv', 'w') as f:
        for i, cs in enumerate(css):
            for c in cs.flatten():
                print(i, c, sep=',', file=f)
    with open('predictions.csv', 'w') as f:
        for i, (ps, test) in enumerate(zip(pss, test_sets)):
            for p, t in zip(ps, test):
                print(i, t, p, y[t], sep=',', file=f)
    rocs = [sorted(roc(ps, y[test])) for ps, test in zip(pss, test_sets)]
    with open('rocs.csv', 'w') as f:
        for i, roc in enumerate(rocs):
            for fpr, tpr, score in roc:
                print(i, fpr, tpr, score, sep=',', file=f)
    with open('avg.csv', 'w') as f:
        for tpr, fpr, _ in sorted(hull(list(avg_roc(rocs)))):
            print(tpr, fpr, sep=',', file=f)
