#!/bin/bash
base=/broad/compbio/aksarkar/projects/impute
diseases=(BD CAD CD HT RA T1D T2D)
cd /broad/compbio/lward/GWAS/ebi_ega_downloads/WTCCC_7_Diseases
for d in ${diseases[@]}
do
    cd $d/Affymetrix500K/Chiamo/OXSTATS_FORMAT
    python $HOME/code/wtccc1/impute/genotypes.py $base/affy/probes.bed $base/genotypes/$d/ *.gz
done
