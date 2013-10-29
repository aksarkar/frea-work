#!/bin/bash
set -e
base=${1?"missing base input dir"}
chain=${2?"missing chain file"}
chr=${3?"missing chromosome"}
zcat $1/*_$(printf '%02d' $chr)_*.gz | awk -vOFS='\t' -vc=$chr '{print "chr"c, $3 - 1, $3, $1}' >hg17.$chr
~lward/bin/liftOver/liftOver hg17.$chr $chain hg19.$chr unmapped.$chr
zcat /broad/compbio/aksarkar/projects/impute2/1kg/ALL.chr${chr}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing.legend.gz | \
    awk -vOFS='\t' -vc=$chr '$5 == "SNP" {print "chr"c, $2 - 1, $2, $1}' | \
    bedtools intersect -a stdin -b hg19.$chr -wb | \
    awk '{print $8, $4, $3}' | \
    sort -k1
rm {hg17,hg19,unmapped}.$chr
