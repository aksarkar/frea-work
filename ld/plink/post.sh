#!/bin/bash
set -e
export LC_ALL=C
join -12 -24 -t$'\t' -o '2.1 2.2 2.3 1.1 2.5' <(zcat $1/*.gz | grep -v CHR | awk -vOFS='\t' '{split($6, snps, "|"); for (i in snps) {print $1":"$2"-"$3, snps[i]}}' | sort -k2) <(zcat $2 | sort -k4) | bedtools sort | gzip
