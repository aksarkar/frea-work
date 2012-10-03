#!/bin/bash
#BSUB -J maf[1-7]
#BSUB -R rusage[mem=6]
#BSUB -o log
#BSUB -q compbio-week
in=$(find /broad/compbio/aksarkar/wtccc1/imputed -type f | sort | sed -n "$LSB_JOBINDEX p")
out=$(basename $in | sed "s/.bed.gz//")
bedtools intersect -a $in -b /broad/hptmp/aksarkar/1kg-maf/ceu-maf.bed.gz -wb | awk -vOFS='\t' '{print $1, $2, $3, $9, $10}' >$out.bed
