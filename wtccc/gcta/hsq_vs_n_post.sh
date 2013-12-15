#!/bin/bash
set -e
pheno=${1?"missing phenotype"}
grep -H "V(G)/Vp_L" hsq/*.hsq | sed "s/:/ /" | awk -vp=$pheno '{sub(/.*\//, "", $1); split($1, a, "."); print a[1], a[2], $3, $4, p}' | sort -k2g >out/$pheno.top.subsets.txt
