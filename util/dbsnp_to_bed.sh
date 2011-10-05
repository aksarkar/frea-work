#!/bin/bash
zcat $1 | egrep "(^rs[[:digit:]]+|GRCh37)" >temp
sed -nre "s/^(rs[[:digit:]]+).*/\1/p" temp >rsids
sed -nre 's/.*chr=([[:digit:]]+) \| chr-pos=([[:digit:]]+).*/chr\1 \2/p' temp | \
    awk -vOFS='\t' '{print $1, $2, $2 + 1}' | \
    paste - rsids
rm temp rsids
