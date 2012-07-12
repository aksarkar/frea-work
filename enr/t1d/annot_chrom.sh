#!/bin/bash
#BSUB -J annot_chrom[1-9]
#BSUB -R rusage[argon_io=3]
#BSUB -o annot_chrom.log
#BSUB -q compbio-week
set -e
celltype=$(sed -n "$LSB_JOBINDEX p" /broad/compbio/aksarkar/dbgap/meta/celltypes.txt | tr [A-Z] [a-z])
zcat /broad/compbio/aksarkar/chromhmm/$celltype.bed.gz | \
    grep $PATTERN | \
    bedtools intersect -a $INFILE -b stdin -c | \
    cut -f6 >$celltype.annot
