#!/bin/bash
export LC_ALL=C
echo "SNP $(awk 'NR > 1 {print $1; print $1}' ../../../data/1kg/ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample | tr '\n' ' ')"
sort -k1 $1 | join - $2 | sort -k2g | cut -d' ' -f1,5-
