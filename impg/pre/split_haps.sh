#!/bin/bash
export LC_ALL=C
echo "SNP $(awk 'NR > 1 {print $1; print $1}' /broad/compbio/aksarkar/projects/impute2/1kg/ALL_1000G_phase1integrated_v3.sample | tr '\n' ' ')"
sort -k1 $1 | join - $2 | sort -k2g | cut -d' ' -f1,5-
