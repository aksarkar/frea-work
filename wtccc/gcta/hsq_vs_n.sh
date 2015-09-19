#!/bin/bash
#
# Estimate h_g^2 for increasing subsets of associated array SNPs in a
# cross-validation setting
#
# Author: Abhishek Sarkar <aksarkar@mit.edu>
set -u
set -e

PREFIX=/broad/hptmp/aksarkar/haplotypes
QC=$PREFIX/T1D-clean
CUTOFFS=(1000 2000 3000 4000 5000 7500 10000 20000 30000 40000 50000 100000 200000)

PCGC="java -jar $HOME/.local/src/PCGCRegression/java-pcgc/PCGCRegression.jar"
PRINT='parallel -j1 echo :::'
export LC_ALL=C

pheno=${1?"missing phenotype"}
p=$(basename $pheno .pheno)
fold=${2?"missing fold"}
MAKE_GRM="gcta64-1.24 --thread-num ${NSLOTS:-1} --bfile $QC --keep $p.$fold.keep --make-grm"
log=$p.$fold.log
if [[ $pheno = "random" ]]
then
    pheno="$PREFIX/T1D.pheno"
    awk '{print $1, $2}' $pheno >$p.$fold.test
    cut -f2 $QC.bim | python $HOME/code/wtccc/gcta/random_order.py $fold >$p.$fold.txt
else
    python $HOME/code/wtccc/gcta/testset.py $fold <$pheno >$p.$fold.test
    plink --memory 4096 --bfile $QC --remove $p.$fold.test --pheno $pheno --1 --logistic --out $p.$fold &>/dev/null
    sort -k9g $p.$fold.assoc.logistic | awk 'NR > 1 {print $2}' >$p.$fold.txt
fi
comm -12 <(sed "s/\t/ /" $PREFIX/T1D.grm.id | sort) <(sort $p.$fold.test) >$p.$fold.keep

$PRINT $p.$fold{,.c} >$p.$fold.grms
$PRINT "$PREFIX/T1D" "$PREFIX/T1D" >$p.$fold.subtractvecs
$PRINT "grmlist=$p.$fold.grms" "subtractvecs=$p.$fold.subtractvecs" "phenos=$pheno" "covars=$PREFIX/T1D.eigenvec" "prevalence=0.005" >$p.$fold.params
for n in ${CUTOFFS[@]}
do
    head -n$n $p.$fold.txt >$p.$fold.head
    $MAKE_GRM --extract $p.$fold.head --out $p.$fold &>$log
    $MAKE_GRM --exclude $p.$fold.head --out $p.$fold.c &>$log
    $PCGC $p.$fold.params >$p.$fold.$n.hsq
done
rm -f $p.$fold.{assoc.logistic,c.*,head,grm*,keep,params,subtractvecs,test,txt}
