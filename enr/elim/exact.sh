#!/bin/bash
#BSUB -J exact
#BSUB -q compbio-week
#BSUB -o out.%J

echo -n "@$ELIM,$CTYPE,$ENH,"
out=$(mktemp -p /broad/shptmp/aksarkar)
$HOME/code/enr/elim/filter $CTYPE $ENH $GWAS $ELIM | cut -f6 >$out
echo $(wc -l <$out) | cat - $out | $HOME/code/enr/exact/ranksum 
rm $out
