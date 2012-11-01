#!/bin/bash
for f in $(find -name "*.out")
do
    pheno=$(echo $f | sed -e 's#./##' -e 's/.bed.out//')
    awk -v p=$pheno '{print $2, $3, p}' $f
done
