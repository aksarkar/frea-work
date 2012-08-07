#!/bin/bash
set -e
ld=/broad/compbio/aksarkar/ld/1kg
bedtools intersect -a stdin -b $ld/haploreg.bed -wb | \
    awk '{print $9, $5}' | \
    python $HOME/code/ld/1kg/expand.py $ld/ceu-index.txt $ld/CEU.txt 0.8 | \
    sort -k1 | \
    join - /broad/compbio/aksarkar/ld/1kg/haploreg.bed -24 | \
    awk -vOFS='\t' '{print $3, $4, $5, $1, $2}'
