#!/bin/bash
#$ -S /bin/bash
#$ -V
#$ -cwd
#$ -j y
#$ -r y
#$ -sync y
#$ -terse
#$ -tc 2000
# Map tasks in parallel to UGER using bounded number of workers
#
# Usage: qsub -t 1-N ... pmap.sh JOBLIST
#
# Reads list of tasks from JOBLIST. Tasks are assigned (divided equally) based
# on job array specification.
#
# Author: Abhishek Sarkar <aksarkar@mit.edu>
set -e
set -u
# UGER munges LD_LIBRARY_PATH and doesn't pass it through even with -V, so fix
# it up here
export LD_LIBRARY_PATH=$LIBRARY_PATH
awk -vn=$SGE_TASK_LAST -vi=$SGE_TASK_ID 'NR % n == i - 1' $1 | parallel --tmpdir /broad/hptmp/aksarkar/tmp --joblog $JOB_NAME.$JOB_ID.$SGE_TASK_ID.joblog -j1 --halt now,fail,1
