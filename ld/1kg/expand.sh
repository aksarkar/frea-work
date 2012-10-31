#!/bin/bash
# Expand SNPs according to LD
# Usage: expand.sh THRESH

# Reads 0-based BED file of SNPs on stdin, writes 0-based BED file on stdout.
# SNPs in LD are identified by the representative they are in LD with.

# Author: Abhishek Sarkar <aksarkar@mit.edu>
set -e
ld=/broad/compbio/aksarkar/ld/1kg
bedtools intersect -a stdin -b $ld/haploreg.bed -wb | \
    awk '{print $9, $5}' | \
    python $HOME/code/ld/1kg/expand.py $ld/ceu-index.txt $ld/CEU.txt $1 | \
    sort -k1 | \
    join - /broad/compbio/aksarkar/ld/1kg/haploreg.bed -24 | \
    sort -k3 | \
    awk -vOFS='\t' '{print $5, $6, $7, $3, $2, $4}'
