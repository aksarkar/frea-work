#!/bin/bash
# Lift over Affy500K genotypes to hg19
# Author: Abhishek Sarkar

# Get the needed Affy500K annotations
[[ ! -f Affy500K_STY_annot.txt.gz ]] && wget ftp://ftp.sanger.ac.uk/pub/WTCCC/misc/Affy500K_STY_annot.txt.gz
[[ ! -f Affy500K_NSP_annot.txt.gz ]] && wget ftp://ftp.sanger.ac.uk/pub/WTCCC/misc/Affy500K_NSP_annot.txt.gz

t=$(mktemp -p /broad/hptmp/aksarkar)
zcat Affy500K_*.txt.gz | grep ^S | awk -vOFS='\t' '{print "chr"$3, $4, $4+1, $1, $6$7, $5}' >$t
~lward/bin/liftOver/liftOver $t ~lward/bin/liftOver/hg17ToHg19.over.chain affy.bed unmapped.log 2>/dev/null
sort -k4 affy.bed >$t
mv $t affy.bed

diseases=(BD CAD CD HT RA T1D T2D)
for d in diseases
do
    for f in $(find /broad/compbio/lward/GWAS/ebi_ega_downloads/WTCCC_7_Diseases/$d/Affymetrix500K/Chiamo/OXSTATS_FORMAT -name "*chiamo.gz")
    do
        zcat $f | sort -k1 | join affy.bed - | cut -f 
