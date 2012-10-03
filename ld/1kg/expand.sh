#!/bin/bash
#BSUB -J expand-wtccc1[1-7]
#BSUB -R rusage[mem=6]
#BSUB -o log
#BSUB -q compbio-week
set -e
in=$(find /broad/compbio/aksarkar/wtccc1/pruned/ -type f | sort | sed -n "$LSB_JOBINDEX p")
out=$(basename $in | sed "s/.gz//")
ld=/broad/compbio/aksarkar/ld/1kg
bedtools intersect -a $in -b $ld/haploreg.bed -wb | \
    awk '{print $9, $5}' | \
    python $HOME/code/ld/1kg/expand.py $ld/ceu-index.txt $ld/CEU.txt $THRESH | \
    sort -k1 | \
    join - /broad/compbio/aksarkar/ld/1kg/haploreg.bed -24 | \
    sort -k3 | \
    awk -vOFS='\t' '{print $5, $6, $7, $3, $2, $4}' >$out
