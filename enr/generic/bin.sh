#!/bin/bash
#BSUB -J wtccc1-preprocess-rrplots[1-1960]
#BSUB -R rusage[argon_io=2]
#BSUB -o wtccc1-preprocess-rrplots.log
#BSUB -q compbio-week
set -e
args=$(sed -n "$LSB_JOBINDEX p" $JOBLIST)
markers=$(echo $args | awk '{print $1}')
features=$(echo $args | awk '{print $2}')
phenotype=$(basename $markers | sed -r 's/.bed(.gz)?//')
f=$(echo $features | sed 's#.*features/##' | awk -vFS=/ '{print $1}')
c=$(echo $features | sed 's#.*features/##' | awk -vFS=/ '{print $2}' | \
    sed -r 's/.bed(.gz)?//')
bedtools intersect -a $markers -b /broad/compbio/aksarkar/wtccc1/exclude.bed.gz -v | \
    bedtools intersect -a stdin -b $features -c | \
    sort -k5g | \
    cut -f6 | \
    python $HOME/code/enr/generic/bin.py $phenotype $c $f 1000 >/broad/hptmp/aksarkar/$(printf "%04d" $LSB_JOBINDEX).out
