#!/bin/bash
set -e
plink --bfile /broad/compbio/aksarkar/projects/gwas/wtccc1/data/genotypes/plink/all --chr $1 --indep-pairwise 300 50 0.2 --out $1
plink --bfile /broad/compbio/aksarkar/projects/gwas/wtccc1/data/genotypes/plink/all --show-tags $1.prune.in --list-all --tag-r2 0.2
