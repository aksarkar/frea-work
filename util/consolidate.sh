#!/bin/bash
# Convert to 0-based and consolidate enhancer dips
# Author: Abhishek Sarkar <aksarkar@mit.edu>
awk -vOFS='\t' '{print $1, $2 - 1, $3}' <$1 | mergeBed -d 1 | bzip2 --best >$1.bz2
