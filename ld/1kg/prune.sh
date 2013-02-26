#!/bin/bash
# Prune SNPs by picking representatives for LD blocks
# Usage: prune.sh THRESH

# Reads 0-based BED file of SNPs on stdin, writes 0-based BED file on stdout

# Author: Abhishek Sarkar <aksarkar@mit.edu>

set -e
ld=/broad/compbio/aksarkar/ld/1kg
bedtools intersect -a stdin -b $ld/haploreg_b137_chromsweep.bed -sorted -wb | \
    awk '{print $9, $5}' | \
    python $HOME/code/ld/1kg/prune.py $ld/ceu-index.txt $ld/CEU.txt $1 | \
    sort -k1 | \
    join - /broad/compbio/aksarkar/ld/1kg/haploreg.bed -24 | \
    awk -vOFS='\t' '{print $3, $4, $5, $1, $2}'
