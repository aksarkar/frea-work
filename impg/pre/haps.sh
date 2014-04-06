#!/bin/bash
# Usage: haps.sh CHROMOSOME
export LC_ALL=C
prefix=../../../data/1kg/ALL.chr${1}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nosing
paste -d' ' <(zcat $prefix.legend.gz | awk 'NR > 1 {print $1}') <(zcat $prefix.haplotypes.gz | tr 01 12 | tr -d ' ') | sort -k1
