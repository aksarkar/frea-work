#!/bin/bash
# Count how many GWAS SNPs fall in enhancer dips
# Author: Abhishek Sarkar <aksarkar@mit.edu>
echo -n $1
bzcat $1 | \
    intersectBed -a $HOME/hp/t1d/data/hg19/filtered-expanded-0.8.bed -b stdin -c | \
    python $HOME/code/ld/group.py union | grep "1$" | wc -l
