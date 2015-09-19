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
awk -vn=$SGE_TASK_LAST -vi=$SGE_TASK_ID 'NR % n == i - 1' $1 | parallel -j1 --halt 1 -t
