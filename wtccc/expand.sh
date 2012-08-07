#!/bin/bash
#BSUB -J expand-wtccc[1-7]
#BSUB -R rusage[mem=5,argon_io=3]
#BSUB -o expand-wtccc.log
#BSUB -q compbio-week
set -e
in=$(find /broad/compbio/aksarkar/wtccc1/markers -type f | \
    sort | \
    sed -n "$LSB_JOBINDEX p")
out=$(basename $in)
zcat $in | \
    $HOME/code/ld/1kg/expand.sh | \
    gzip >$out
