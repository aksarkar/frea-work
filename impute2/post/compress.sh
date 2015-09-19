#!/bin/bash
#$ -pe smp 8
set -e
set -u
cohort=${1?"missing cohort"}
parallel --halt 1 "cat ${cohort}_{}_*.gen | gzip >${cohort}_{}.gen.gz && mv ${cohort}_{}.gen.gz /broad/hptmp/aksarkar/out" ::: $(seq -w 1 22)
parallel --halt 1 -X rm ::: ${cohort}_* impute2.*
