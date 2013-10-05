#!/bin/bash
set -e
plink --tped $1 --tfam $2 --chr $5 --update-map $3/pos.txt --make-bed --out $4.pos
plink --bfile $4.pos --update-map $3/names.txt --update-name --extract $3/rsids.txt --mind .01 --geno .01 --maf .01 --hwe .05 --make-bed --out $4
