#!/bin/bash
#BSUB -J dbgap-dnase[1-4215]
#BSUB -R rusage[argon_io=3]
#BSUB -o dbgap-dnase.log
#BSUB -q compbio-week
set -e
dbgap=/broad/compbio/aksarkar/dbgap
markers=$(find $dbgap/analyses -type f | \
    sort | \
    sed -n "$LSB_JOBINDEX p")
expanded=$dbgap/expanded/hapmap/$(printf "%04d" $LSB_JOBINDEX).bed.gz
annotations=/broad/compbio/aksarkar/annotations
for i in $(find $annotations/dnase/ -type f)
do
    echo -n ">>> $markers $(echo $i | sed -r 's#.*/(.*).bed.gz#\1#') "
    bedtools intersect -a $expanded -b $i -c | \
        python3 $HOME/code/ld/group.py union | \
        sort -g -k5 | \
        cut -f5,7 | \
        python $HOME/code/test/ks.py
done
