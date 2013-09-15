#!/bin/bash
export LC_ALL=C
awk 'NR > 1' ${1?"missing"} | sort -k2 | awk '{print $1, $4, $14 / $15 >$2".tmp"}'
seq 1 22 | parallel 'LC_ALL=C sort -k2 ../maps/{}.map | LC_ALL=C join {}.tmp - -12 -22 -o "1.1 1.2 2.3 2.4 1.3" | sort -k2g | cat <(echo "snp pos ref alt z") - >{}.txt'
rm -f *.tmp
