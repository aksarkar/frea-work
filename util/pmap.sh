#!/bin/bash
#BSUB -g /compbio/short
#BSUB -o log
#BSUB -q compbio-week

# Map tasks in parallel to LSF using bounded number of workers 

# Reads list of tasks from file specified in environment variable JOBLIST.
# Workers are assigned based on LSF job array specification.

# Author: Abhishek Sarkar <aksarkar@mit.edu>

set -e
n=$(wc -l <$JOBLIST)
LSB_JOBINDEX_START=$(echo $LSB_JOBNAME | sed -e 's/.*\[//' -e 's/-.*//')
m=$((($LSB_JOBINDEX_END - $LSB_JOBINDEX_START) / $LSB_JOBINDEX_STEP + 1))
if (($m > $n))
then
    exit 1
fi
ntasks=$(($n / $m))
rem=$(($n % $m))
if (($LSB_JOBINDEX - 1 < $rem))
then
    ((ntasks++))
fi
start=$(($ntasks * ($LSB_JOBINDEX - 1) + 1))
if (($LSB_JOBINDEX > $rem))
then
    start=$(($start + $rem))
fi
end=$(($LSB_JOBINDEX == $LSB_JOBINDEX_END ? $n + 1 : $start + $ntasks))
for ((i = $start; i < $end; i++))
do
    echo "$(sed -n ${i}p $JOBLIST)"
    eval "$(sed -n ${i}p $JOBLIST)"
done
