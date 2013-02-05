#!/bin/bash
#
# Compute the proportion of SNPs in LD with the specified regions at each LD
# block size
#
# Usage: expectation.sh REGIONS THRESH
#   REGIONS - 0-based BED
#   THRESH - R^2 threshold
#
# Author: Abhishek Sarkar <aksarkar@mit.edu>
ld=/broad/compbio/aksarkar/ld/1kg
bedtools intersect -a $ld/haploreg.bed -b $1 -c | \
    python $HOME/code/enr/generic/expectation/expectation.py $2 $ld/CEU.txt
