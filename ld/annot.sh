#!/bin/bash

# Annotate GWAS markers by expanding according to LD and intersecting with
# regions

# Author: Abhishek Sarkar <aksarkar@mit.edu>

# Usage: annot.sh THRESH REGIONS FN
# THRESH - minimum R^2 to consider
# REGIONS - BED file of regions
# FN - function to combine annotations. One of 'state' or 'union'

# Reads GWAS BED file on stdin. Writes annotated BED file to stdout

python $HOME/code/ld/lookup.py $1 | \
    intersectBed -a stdin -b $2 -c | \
    python $HOME/code/ld/group.py $3
