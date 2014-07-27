"""Write motif co-occurence network for GraphViz

Usage: python network.py OUTFILE THRESH

THRESH is minimum proportion of co-occurences to individual occurences (used to
prune directed edges).

Expects space-separated (chr, start, end, motif) on stdin. Position information
is used to determine co-occurence per region.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
from __future__ import print_function

import collections
import operator
import itertools
import sys

import pygraphviz

thresh = float(sys.argv[2])

pairs = collections.Counter()
data = (line.split() for line in sys.stdin)
for k, g in itertools.groupby(data, key=operator.itemgetter(0,1,2)):
    motifs = [x[3] for x in g]
    pairs.update((x, y) for x in motifs for y in motifs)

graph = pygraphviz.AGraph(strict=False, directed=True)
graph.graph_attr['size'] = '4,4'
graph.graph_attr['start'] = 0
graph.graph_attr['maxiter'] = 50000
graph.graph_attr['splines'] = 'curved'
graph.node_attr['fillcolor'] = '#ffffffaa'
graph.node_attr['fontname'] = 'Nimbus Sans L'
graph.node_attr['shape'] = 'rectangle'
graph.node_attr['style'] = 'filled',
graph.edge_attr['color'] = '#00000020'
graph.edge_attr['arrowSize'] = .5
for i, j in pairs:
    if i != j:
        n = pairs[(i, j)]
        if n >= thresh * pairs[(i, i)]:
            graph.add_edge(i, j)
graph.write(path=sys.argv[1])
