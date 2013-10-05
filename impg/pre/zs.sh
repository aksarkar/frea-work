#!/bin/bash
export LC_ALL=C
scores=${1?"missing scores"}
chrom=${2?"missing chromosome"}
echo "snp pos ref alt z"
awk -vc=$chrom '$2 == c' $scores | \
    sort -k2 | \
    uniq | \
    awk '{print $1, $4, $14 / $15}' | \
    join - <(sort -k2 ../maps/$chrom.map) -12 -22 -o "1.1 1.2 2.3 2.4 1.3" | \
    awk '$3 ~ /^.$/ && $4 ~ /^.$/' | \
    sort -k2g
