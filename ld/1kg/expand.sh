#!/bin/bash
# Expand SNPs according to LD
# Usage: expand.sh THRESH

# Reads 0-based BED file of SNPs on stdin, writes 0-based BED file on stdout.
# SNPs in LD are identified by the representative they are in LD with.

# Author: Abhishek Sarkar <aksarkar@mit.edu>
set -e
export LC_ALL=C
ld=/broad/compbio/aksarkar/ld/1kg
bedtools intersect -a $ld/haploreg_b137_chromsweep.bed.gz -b stdin -sorted | \
    cut -f4 | \
    sort | \
    join - $ld/normalized.txt | \
    awk -v thresh=$1 '$3 > thresh {print $2}' | \
    sort | \
    uniq | \
    join - $ld/haploreg.bed -24 -o '2.1 2.2 2.3 2.4' | \
    tr ' ' '\t' | \
    bedtools sort
