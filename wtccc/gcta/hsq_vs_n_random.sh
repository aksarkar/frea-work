#!/bin/bash
set -e
pheno=${1?"missing phenotype"}
fold=${2?"missing fold"}
n="5000 10000 20000 30000 40000 50000 100000 271428"
qc=/broad/compbio/aksarkar/projects/gwas/wtccc1/results/gcta/2013-12-06/qc
cut -f2 $qc.bim | shuf >raw/$fold.txt
parallel -j1 'gcta64 --thread-num 16 --bfile '$qc' --extract <(head -n{2} raw/{1}.txt) --make-grm-bin --out raw/{1}.{2} &>raw/{1}.{2}.grm.log' ::: $fold ::: $n
parallel -j1 'gcta64 --thread-num 16 --grm-bin raw/{1}.{2} --grm-cutoff .05 --reml --pheno '$pheno' --prevalence .005 --out raw/{1}.{2} &>raw/{1}.{2}.hsq.log' ::: $fold ::: $n
