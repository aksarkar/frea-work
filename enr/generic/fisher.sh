#!/bin/bash
#BSUB -J wtccc1-fisher[1-1960]
#BSUB -R rusage[argon_io=1]
#BSUB -o wtccc1-fisher.log
#BSUB -q compbio-week
set -e
. $HOME/py32/bin/activate
args=$(sed -n "$LSB_JOBINDEX p" $JOBLIST)
markers=$(echo $args | awk '{print $1}')
features=$(echo $args | awk '{print $2}')
bedtools intersect -a $markers -b /broad/compbio/aksarkar/wtccc1/exclude.bed.gz -v | \
    bedtools intersect -a stdin -b /broad/compbio/aksarkar/annotations/ccds/gene.bed.gz -v | \
    bedtools intersect -a stdin -b $features -c | \
    cut -f5,6 | \
    python $HOME/code/test/fisher.py $(basename $markers) $(echo $features | sed -r -e "s#[/a-z]*features/##" -e "s/.bed(.gz)?//")
