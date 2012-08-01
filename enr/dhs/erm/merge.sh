#!/bin/bash
#BSUB -J merge-erm
#BSUB -R rusage[mem=6,argon_io=3]
#BSUB -o merge-erm-dhs.log
#BSUB -q compbio-week
set -e
t=$(mktemp -p /broad/hptmp/aksarkar)
for d in $(find /broad/compbio/eaton/erm_data/Complete-Epigenomes -name "Chromatin_Accessibility" | \
    sed -r 's#/(Donor|Chromatin).*##' | \
    sort | \
    uniq)
do
    out=$(echo $d | sed -r 's#.*/##').bed
    for f in $(find $d -wholename "*/Chromatin_Accessibility/*.bed.gz")
    do
        if [[ ! -f $out ]]
        then
            bedtools merge -i $f >$out
        else
            bedtools merge -i $f | \
                bedtools intersect -a $out -b stdin >$t
            mv $t $out
        fi
    done
    gzip $out
done
