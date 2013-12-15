#!/bin/bash
export LC_ALL=C
markers=${1?"$0: missing markers"}
features=${2?"$0: missing features"}
binsize=${3-1000}
phenotype=$(basename $markers .bed.gz)
f=$(echo $features | sed "s#.*features/##" | cut -d/ -f1)${mask+$mod$(basename $mask | sed "s/.bed.*//")}
c=$(echo $features | sed "s#.*features/##" | cut -d/ -f2- | \
    sed -r "s/.bed.*//")
join <(zcat $markers | sort -k1) <(zcat $features) -o '2.2 2.3 1.2 2.2 2.3 1.2' | \
    sort -k1 | \
    awk -f $HOME/code/ld/union.awk | \
    sort -k1gr | \
    cut -f2 | \
    python $HOME/code/enr/generic/bin/bin.py $phenotype $f $c $binsize 1
