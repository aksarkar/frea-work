#!/bin/bash
for f in $(find /broad/compbio/aksarkar/wtccc1/imputed -name "*.bed.gz")
do
    echo "zcat $f | $HOME/code/ld/1kg/prune.sh .2 >$(basename $f | sed s/.bed.gz//).out"
done
