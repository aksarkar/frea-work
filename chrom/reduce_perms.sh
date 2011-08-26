#!/bin/bash
#BSUB -J reduce_perms
#BSUB -o /dev/null
#BSUB -q compbio-week
work="/broad/shptmp/aksarkar/"
echo "cell_type,state,p" > $PT_OUT
cat $PT_WORK/?? | paste -d, "$PT_WORK/joblist" - >> $PT_OUT
