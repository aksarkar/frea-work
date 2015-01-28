"""Extract known protein-coding genes from GENCODE

Usage: python parse_gtf.py

Expects GTF lines on stdin. Writes BED file with ENSEMBL IDs and gene names on
stdout.

Author: Abhishek Sarkar <aksarkar@mit.edu>

"""
import re
import signal
import sys

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

ensembl_id = re.compile('ENSG[0-9]+')

filter_comment = (line for line in sys.stdin if line[0] != '#')
filter_track = (line for line in filter_comment if not line.startswith('track'))
data = (line.split('\t') for line in filter_track)
for row in data:
    chromosome, source, feature, start, end, score, strand, frame, attribute = row
    attrs = dict(tuple(x.strip().split(' ')) for x in attribute.strip().split(';')[:-1])
    if (source == 'HAVANA' and feature == 'gene' and 
        attrs['gene_type'] == '"protein_coding"' and
        attrs['gene_status'] == '"KNOWN"'):
        gene_id = ensembl_id.search(attrs['gene_id'])
        if gene_id is not None:
            print(chromosome, start, end, gene_id.group(), attrs['gene_name'][1:-1], sep='\t')
