#!/bin/bash
set -e
pheno=${1?"missing phenotype"}
fold=${2?"missing fold"}
n="100 500 1000 5000 10000 50000 100000 271428"
qc=/broad/compbio/aksarkar/projects/gwas/wtccc1/results/gcta/2013-12-06/qc
plink --bfile $qc --remove test-sets/$fold.txt --pheno $pheno --logistic --out rankings/$fold &>/dev/null
sort -k9g rankings/$fold.assoc.logistic | awk 'NR > 1 {print $2}' >rankings/$fold.txt
parallel -j1 'gcta64 --thread-num 16 --bfile '$qc' --keep test-sets/{1}.txt --extract <(head -n{2} rankings/{1}.txt) --maf .01 --make-grm-bin --out grms/{1}.{2} &>grms/{1}.{2}.log' ::: $fold ::: $n
parallel -j1 'gcta64 --thread-num 16 --grm-bin grms/{1}.{2} --grm-cutoff .05 --reml --pheno '$pheno' --prevalence .005 --out hsq/{1}.{2} &>hsq/{1}.{2}.log' ::: $fold ::: $n
