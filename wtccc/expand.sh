#!/bin/bash
set -e
set -u
join -13 <(sort -k3 $1.$2.names) <(sort -k1 $2.expanded) -o "2.2 2.1 1.2" | \
    sort -k1,1 -k3,3gr | \
    awk 'NR==1{k=$1}$1!=k{print;k=$1}END{print}' | \
    join - $2.positions -o "2.2 2.3 1.2 1.3" | \
    awk -vOFS="\t" '{print $1, $2-1, $2, $3, $4}'
