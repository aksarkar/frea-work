#!/bin/bash
#
# Estimate h_g^2 for increasing subsets of associated array SNPs in a
# cross-validation setting
#
# Author: Abhishek Sarkar <aksarkar@mit.edu>
set -e
set -u

COHORT=("train" "test")
CUTOFFS=(1000 2000 3000 4000 5000 7500 10000 20000 30000 40000 50000 100000 200000)

CODE=$HOME/code/wtccc/gcta
PCGC="java -jar $HOME/.local/src/PCGCRegression/java-pcgc/PCGCRegression.jar"
PRINT='parallel -j1 echo :::'
export LC_ALL=C

PREFIX=${1?"missing prefix"}
pheno=$PREFIX.pheno
p=$(basename $pheno .pheno)
seed=${2:-0}
MAKE_GRM="gcta64-1.24 --thread-num ${NSLOTS:-1} --bfile $PREFIX-clean --make-grm"
log=$p.$seed.log

python -c "import sys; assert sys.version_info[0] == 3"
python $CODE/testset.py $seed <$pheno >$p.keep
comm -12 <(sed "s/\t/ /" $PREFIX.grm.id | sort) <(sort $p.keep) >$p.train
comm -23 <(sed "s/\t/ /" $PREFIX.grm.id | sort) <(sort $p.keep) >$p.test

$PRINT $p{,.c} >$p.grms
$PRINT "$PREFIX" "$PREFIX" >$p.subtractvecs
$PRINT "grmlist=$p.grms" "subtractvecs=$p.subtractvecs" "phenos=$pheno" "covars=$PREFIX.eigenvec" "prevalence=0.005" >$p.params
for c in ${COHORT[@]}
do
    if [[ seed == 0 ]]
    then
        plink --memory 4096 --bfile $PREFIX-clean --remove $p.$c --pheno $pheno --1 --logistic --out $p &>/dev/null
        sort -k9g $p.assoc.logistic | awk 'NR > 1 {print $2}' >$p.txt
    else
        cut -f2 $PREFIX-clean.bim | python $CODE/random_order.py $seed >$p.txt
    fi
    for n in ${CUTOFFS[@]}
    do
        head -n$n $p.txt >$p.head
        $MAKE_GRM --keep $p.$c --extract $p.head --out $p &>>$log
        $MAKE_GRM --keep $p.$c --exclude $p.head --out $p.c &>>$log
        if [[ seed == 0 ]]
        then
            hsq=$p.$n.$c.hsq
        else
            hsq=$p.$n.$c-$seed.hsq
        fi            
        $PCGC $p.params >$hsq 2>>$log
    done
done
rm -f $p.{assoc.logistic,c.*,head,grm*,keep,params,subtractvecs,train,test,txt}
