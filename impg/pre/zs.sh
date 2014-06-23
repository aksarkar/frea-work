#!/bin/bash
export LC_ALL=C
scores=${1?"missing scores"}
chrom=${2?"missing chromosome"}
echo "snp pos ref alt z"
awk -vc=$chrom '$2 == c' $scores | \
    sort -k2 | \
    uniq | \
    awk '{print $1, $4, "A", "C", $14 / $15}' | \
    sort -k2g
