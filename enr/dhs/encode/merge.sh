#!/bin/bash
# Intersect DNAse replicates
# Author: Abhishek Sarkar <aksarkar@mit.edu>
set -e
t=$(mktemp -p /broad/hptmp/aksarkar)
export DNASE=/broad/compbio/rca/to_sort/compbio-t_data/encode/UwDnase/
sed -e "s/; / /g; s/cell=//" $DNASE/files.txt | \
    awk '/narrowPeak/ {print $5, $1}' >$t
awk '{print $1}' $t | uniq | \
    parallel $HOME/code/util/helper.sh $t
rm $t
