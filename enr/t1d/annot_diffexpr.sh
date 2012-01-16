#!/bin/bash
# Usage: annot_diffexpr.sh GWAS

ctypes=(gm12878 h1 hepg2 hmec hsmm huvec k562 nhek nhlf)
enh=(Strong_Enhancer Weak_Enhancer Enhancer)

annot=~pouyak/work/motifs/verts/stats/encode/brad/special/diff-gene-regions.txt

for c in ${ctypes[@]}
do
    awk -vOFS='\t' '($1 == "CDS" || $1 == "3UTR" || $1 == "5UTR") && $7 == "protein_coding" && tolower($8) == "'$c'" {print $2, $3 - 1, $4}' $annot | \
        intersectBed -a $1 -b stdin -c | cut -f6>$c.annot
done
echo "rsid,p,gm12878,h1,hepg2,hmec,hsmm,huvec,k562,nhek,nhlf"
zcat $1 | awk -vOFS=',' '{print $4, $5}' | paste -d, - *.annot | sort -t, -g -k2 | grep -v ",,"
rm *.annot
