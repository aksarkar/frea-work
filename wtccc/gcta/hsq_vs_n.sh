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
prevalence=${2?"missing prevalence"}
mode=${3?"missing mode"}
c=${4?"missing hold-out chromosome"}

if [[ observed == ${mode}* ]]
then
    mode=observed
elif [[ random == ${mode}* ]]
then
    mode=random
else
    exit 1
fi
p=$(basename $PREFIX)
mkdir -p $p.$mode.$c
pushd $p.$mode.$c

python -c "import sys; assert sys.version_info[0] == 3"
if [[ $mode == "observed" ]]
then
   python $CODE/testset.py 0 <$PREFIX.pheno >$p.keep
   comm -12 <(sed "s/\t/ /" $PREFIX.grm.id | sort) <(sort $p.keep) >$p.train
   comm -23 <(sed "s/\t/ /" $PREFIX.grm.id | sort) <(sort $p.keep) >$p.test
else
    parallel -j1 ln -sf $PREFIX-clean.{} $p.{} ::: bed bim fam
    ln -sf $PREFIX.grm.id $p.test
    cut -f2 $PREFIX-clean.bim | python $CODE/random_order.py 0 >$p.txt
fi

$PRINT $p{,.c} >$p.grms
$PRINT "$PREFIX" "$PREFIX" >$p.subtractvecs
$PRINT "grmlist=$p.grms" "subtractvecs=$p.subtractvecs" "phenos=$PREFIX.pheno" "covars=$PREFIX.eigenvec" "prevalence=$prevalence" >$p.params

PLINK="plink --memory 4096 --out $p"
for n in ${CUTOFFS[@]}
do
    $PLINK --bfile $PREFIX-clean --exclude range <(echo "$c 0 100000000 drop") --make-bed &>/dev/null
    if [[ $mode == "observed" ]]
    then
        $PLINK --bfile $p --remove $p.test --pheno $PREFIX.pheno --1 --logistic --covar $PREFIX.eigenvec &>/dev/null
        sort -k9g $p.assoc.logistic | awk 'NR > 1 {print $2}' >$p.txt
    fi
    MAKE_GRM="gcta64-1.24 --thread-num 8 --bfile $p --keep $p.test --make-grm"
    head -n$n $p.txt >$p.head
    $MAKE_GRM --extract $p.head --out $p &>>log
    $MAKE_GRM --exclude $p.head --out $p.c &>>log
    $PCGC $p.params >$p.$c.$n.hsq 2>>log
done
rm -f $p.{assoc.logistic,bed,bim,fam,c.*,head,grm*,keep,params,subtractvecs,train,test,txt}
