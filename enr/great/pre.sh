#!/bin/bash
bedtools intersect -a 30k.bed.gz -b $1 -sorted -wao | \
    awk -vOFS='\t' '$6 != "." {print $6, $7, $8}'
