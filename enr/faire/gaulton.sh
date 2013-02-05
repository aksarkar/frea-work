#!/bin/bash
th=(liberal moderate stringent)
for t in ${th[@]}
do
    a=""
    for f in $(find -name "*$t*")
    do
        if [[ -z $a ]]
        then
            echo -n "cat $f | "
            a=$f
        else
            echo -n "bedtools intersect -a stdin -b $f | "
        fi
    done
    echo "gzip >$t.bed.gz"
done
