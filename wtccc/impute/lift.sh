#!/bin/bash
set -e
c=$(basename $1 | awk -vFS='.' '{sub(/^0/, "", $2); print $2}')
hg17=$(basename $1 .gen.gz).hg17.bed
hg19=$(echo $hg17 | sed "s/hg17/hg19/")
snps=$(basename $1 .gen.gz).txt
[[ ! -f $hg17 ]] && zcat $1 | awk -vOFS='\t' -vc=$c '{print "chr"c, $3 - 1, $3, $3 - 1}' >$hg17
[[ ! -f $hg19 ]] && ~lward/bin/liftOver/liftOver $hg17 ~lward/bin/liftOver/hg17ToHg19.over.chain $hg19 $(basename hg17 .bed).log
[[ ! -f $snps ]] && awk 'NF == 4' $hg19 | \
    bedtools sort | \
    bedtools intersect -a stdin -b /broad/compbio/aksarkar/ld/1kg/haploreg_b137_chromsweep.bed.gz -sorted -wb | \
    awk '{sub(/chr/, "", $1); print $1, $4, $2, $8}' >$snps
zcat $1 | python $HOME/code/wtccc/impute/lift.py $snps >$(basename $1 .gz)
