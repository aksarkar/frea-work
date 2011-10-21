#!/bin/bash
# Intersect DNAse replicates
# Author: Abhishek Sarkar <aksarkar@mit.edu>
t=$(mktemp -p $HOME/tmp)
sed -e "s/; / /g; s/cell=//" /fg/compbio-t/rca/data/encode/UwDnase/files.txt | \
    awk '/narrowPeak/ {print $5, $1}' >$t
awk '{print $1}' $t | uniq | \
    parallel $(cd -P $(dirname $0) && pwd)/merge_dnase_helper.sh $t
rm $t
