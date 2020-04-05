#!/bin/bash
set -e
pheno=${1?"missing phenotype"}
p=$(basename $pheno .txt)
fold=${2?"missing fold"}
PREFIX=/broad/compbio/aksarkar/projects/gwas/wtccc1/EC21
qc=$PREFIX/results/qc/merge/qc.gcta
prune=$PREFIX/results/gcta/h2g-all-samples/prune.txt
n="1000 2000 3000 4000 5000 7500 10000 20000 30000 40000 50000 100000 271428"
if [[ $pheno = "random" ]]
then
    cut -f2 $qc.bim | python $HOME/code/gcta/random_order.py $fold >rankings/$p.$fold.txt
    parallel -j1 "gcta64-1.24 --bfile $qc --remove $prune --extract "'<('"head -n{} rankings/$p.$fold.txt"')'" --make-grm --out grms/$p.$fold.{} &>grms/$p.$fold.{}.log" ::: $n
    pheno="$PREFIX/results/gcta/h2g-all-samples/t1d.txt"
else
    python $HOME/code/gcta/testset.py $fold <$pheno >test-sets/$p.$fold.txt
    $PLINK --bfile $qc --remove test-sets/$p.$fold.txt --pheno $pheno --1 --logistic --out rankings/$p.$fold &>/dev/null
    sort -k9g rankings/$p.$fold.assoc.logistic | awk 'NR > 1 {print $2}' >rankings/$p.$fold.txt
    parallel -j1 "gcta64-1.24 --bfile $qc --keep test-sets/$p.$fold.txt --remove $prune --extract "'<('"head -n{} rankings/$p.$fold.txt"')'" --make-grm --out grms/$p.$fold.{} &>grms/$p.$fold.{}.log" ::: $n
fi
parallel -j1 "gcta64-1.24 --grm grms/$p.$fold.{} --reml --pheno $pheno --prevalence .005 --out hsq/$p.$fold.{} &>hsq/$p.$fold.{}.log" ::: $n
