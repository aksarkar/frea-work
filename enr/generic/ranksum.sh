#!/bin/bash
#BSUB -J wtccc1-ranksum[1-1960]
#BSUB -R rusage[argon_io=1]
#BSUB -o wtccc1-generic.log
#BSUB -q compbio-week
set -e
args=$(sed -n "$LSB_JOBINDEX p" $JOBLIST)
markers=$(echo $args | awk '{print $1}')
features=$(echo $args | awk '{print $2}')
echo -n ">>> $(basename $markers) $features " | \
    sed -r -e "s#[/a-z]*features/##" -e "s/.bed(.gz)?//"
t=$(mktemp -p /broad/hptmp/aksarkar)
bedtools intersect -a $markers -b /broad/compbio/aksarkar/wtccc1/exclude.bed.gz -v | \
    bedtools intersect -a stdin -b $features -c | \
    bedtools sort -i stdin | \
    cut -f5,6 >$t
wc -l <$t | cat - $t | $HOME/code/test/exact/ranksum
rm $t
