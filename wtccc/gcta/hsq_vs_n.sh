#!/bin/bash
set -e
pheno=${1?"missing phenotype"}
p=$(basename $pheno .txt)
fold=${2?"missing fold"}
qc=/broad/compbio/aksarkar/projects/gwas/wtccc1/results/gcta/2013-12-06/qc
prune=/broad/compbio/aksarkar/projects/gwas/wtccc1/results/gcta/2013-12-06/h2g-all-samples/prune.txt
python $HOME/code/wtccc/gcta/testsets.py <$pheno >test-sets/$p.$fold.txt
plink --bfile $qc --remove test-sets/$p.$fold.txt --pheno $pheno --1 --logistic --out rankings/$p.$fold &>/dev/null
sort -k9g rankings/$p.$fold.assoc.logistic | awk 'NR > 1 {print $2}' >rankings/$p.$fold.txt
n="1000 2000 3000 4000 5000 7500 10000 25000 50000 100000 271428"
parallel -j1 "gcta64 --thread-num 16 --bfile $qc --keep test-sets/$p.$fold.txt --remove $prune --extract "'<('"head -n{} rankings/$p.$fold.txt"')'" --make-grm --out grms/$p.$fold.{} &>grms/$p.$fold.{}.log" ::: $n
parallel -j1 "gcta64 --thread-num 16 --grm grms/$p.$fold.{} --reml --pheno $pheno --prevalence .005 --out hsq/$p.$fold.{} &>hsq/$p.$fold.{}.log" ::: $n
