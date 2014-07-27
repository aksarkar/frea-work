#!/bin/bash
set -e
scores=${1?"missing scores"}
zcat $scores | sort -k5gr | awk -v s=$scores 'BEGIN {n = 100} NR == n {print s, $5; n *= 2}'
