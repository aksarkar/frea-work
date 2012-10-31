#!/bin/bash
# Bin binary annotations for RR plots
# Usage: bin.sh MARKERS FEATURES BINSIZE [THRESH]

#   MARKERS - 0-based BED file. If THRESH is given, 4th field must be
#   representative id.

#   FEATURES - 0-based BED file

#   BINSIZE - number of SNPs per bin

#   THRESH - LD threshold to take union over

# Author: Abhishek Sarkar <aksarkar@mit.edu>
set -e
markers=$1
features=$2
phenotype=$(basename $markers | sed -r 's/.bed(.gz)?//')
f=$(echo $features | sed 's#.*features/##' | awk -vFS=/ '{print $1}')
c=$(echo $features | sed 's#.*features/##' | awk -vFS=/ '{print $2}' | \
    sed -r 's/.bed(.gz)?//')
if [[ -z $4 ]]
then
    bedtools intersect -a $1 -b $features -c | \
        sort -k5g | \
        cut -f6 | \
        python $HOME/code/enr/generic/bin/bin.py $phenotype $f $c $3
else
    zcat $markers | \
        awk -vt=$4 '$6 > t' | \
        bedtools intersect -a stdin -b $features -c | \
        python $HOME/code/ld/group.py union | \
        sort -k5g | \
        cut -f7 | \
        python $HOME/code/enr/generic/bin/bin.py $phenotype $f $c $3
fi
