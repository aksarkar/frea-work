#!/bin/bash
#BSUB -J postprocess[1-4215]
#BSUB -R rusage[argon_io=3]
#BSUB -o log
#BSUB -q compbio-week
set -e
i=$(printf "%04d" $LSB_JOBINDEX)
zcat /broad/compbio/aksarkar/dbgap/for-luke/$i.out.gz | \
    sort -k1 | \
    join - /broad/compbio/aksarkar/ld/1kg/haploreg.bed -24 | \
    awk -vOFS='\t' '{print $3, $4, $5, $1, $2}' | \
    gzip >$i.bed.gz
