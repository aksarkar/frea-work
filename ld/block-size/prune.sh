#!/bin/bash
# Prune SNPs by picking representatives for LD blocks
# Usage: prune.sh THRESH

# Reads 0-based BED file of SNPs on stdin, writes 0-based BED file on stdout

# Author: Abhishek Sarkar <aksarkar@mit.edu>

set -e
ld=/broad/compbio/aksarkar/ld/1kg
bedtools intersect -a stdin -b $ld/haploreg.bed -wb | \
    awk '{print $9, $5}' | \
    python $HOME/code/ld/1kg/prune.py $ld/ceu-index.txt $ld/CEU.txt $1
