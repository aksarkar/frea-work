from __future__ import print_function

import csv
import os
import sys

import suds.client

with open(sys.argv[1]) as f:
    ids = ','.join(line.strip() for line in f)
with open(sys.argv[2]) as f:
    bg = ','.join(line.strip() for line in f)
if not ids:
    sys.exit(0)
david = suds.client.Client('http://david.abcc.ncifcrf.gov/webservice/services/DAVIDWebService?wsdl').service
david.authenticate(b'aksarkar@mit.edu')
epigenome, _ = os.path.splitext(os.path.basename(sys.argv[1]))
david.addList(ids, 'ENSEMBL_GENE_ID', epigenome, 0)
david.addList(bg, 'ENSEMBL_GENE_ID', "back", 1)
clusters = david.getTermClusterReport(4, 4, 4, .5, 35)
w = csv.writer(sys.stdout)
for i, cluster in enumerate(clusters):
    for term in cluster.simpleChartRecords:
        w.writerow([i, cluster.score, term.categoryName, term.termName, '{:.3g}'.format(term.ease), '{:.3g}'.format(term.foldEnrichment)])
