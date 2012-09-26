#!/bin/bash
#BSUB -J intersect-wtccc1
#BSUB -R rusage[argon_io=1]
#BSUB -o intersect-wtccc1.log
#BSUB -q compbio-week
set -e
echo "study imputed expanded" >out.txt
for i in {1..7}
do
    imputed=$(find /broad/compbio/aksarkar/wtccc1/imputed -type f | \
        sort | \
        sed -n "$i p")
    expanded=$(find /broad/compbio/aksarkar/wtccc1/expanded -type f | \
        sort | \
        sed -n "$i p")
    study=$(basename $imputed | sed 's/.bed.gz//')
    bedtools intersect -a $imputed -b /broad/compbio/aksarkar/wtccc1/exclude.bed.gz -v | \
        bedtools intersect -a stdin -b $expanded -wb | \
        awk -vs=$study '{print s, $5, $10}' >>out.txt
done
