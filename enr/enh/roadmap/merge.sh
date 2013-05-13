#!/bin/bash
#BSUB -J merge-erm
#BSUB -R rusage[mem=1,argon_io=3]
#BSUB -o merge-erm-enh.log
#BSUB -q compbio-week
set -e
t=$(mktemp -p /broad/hptmp/aksarkar)
for f in $(find -name "*.bed.gz" | \
    sed -r 's#(.Donor.*)?.bed.gz##' | \
    sort | \
    uniq)
do
    out=/broad/hptmp/aksarkar/$(echo $f | sed -r 's#.*/##').bed
    for g in $(find -wholename "$f*")
    do
        if [[ ! -f $out ]]
        then
            bedtools merge -i $g >$out
        else
            bedtools merge -i $g | \
                bedtools intersect -a $out -b stdin -sorted >$t
            mv $t $out
        fi
    done
    # gzip $out
done
