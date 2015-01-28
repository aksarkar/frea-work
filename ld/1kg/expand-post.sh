#!/bin/bash
set -e
set -u
parallel -j1 'join -14 -22 <(zcat chr{}.txt.gz | sort -k4) <(zcat chr{}-pruned.gz | sort -k2)' ::: $(seq 1 22) | awk -vOFS='\t' '$6 == "" {$6 = $1} {print $2, $3, $4, $6, $5}' | bedtools sort
