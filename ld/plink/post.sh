#!/bin/bash
export LC_ALL=C
awk -f $HOME/code/ld/plink/post.awk | \
    sort -k2 | \
    join - t1d.txt -12 -24 -t$'\t' -o '2.1 2.2 2.3 1.1 2.5' | \
    bedtools sort
