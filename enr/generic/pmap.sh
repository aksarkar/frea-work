#!/bin/bash
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
start=$(($n / $m * ($LSB_JOBINDEX - 1) + 1))
end=$(($LSB_JOBINDEX == $LSB_JOBINDEX_END ? $n + 1 : $start + $n / $m))
for ((i = $start; i < $end; i++))
do
    eval $(sed -n "$i p" $JOBLIST)
done
