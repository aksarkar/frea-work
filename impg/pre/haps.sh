#!/bin/bash
# Usage: haps.sh CHROMOSOME
export LC_ALL=C
paste -d' ' <(zcat /broad/compbio/aksarkar/projects/impute2/1kg/ALL_1000G_phase1integrated_v3_chr${1}_impute_macGT1.legend.gz | awk 'NR > 1 {print $1}') <(zcat /broad/compbio/aksarkar/projects/impute2/1kg/ALL_1000G_phase1integrated_v3_chr${1}_impute_macGT1.hap.gz | tr 01 12 | tr -d ' ') | sort -k1
