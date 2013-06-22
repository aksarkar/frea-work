#!/bin/bash
awk -vc=${2:-.01} -vFS='\t' -vOFS='\t' '$1 == "GO Biological Process" && $7 < c && $16 < c && $8 > 1 {print FILENAME, $3, $8}' $1
