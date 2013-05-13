#!/bin/bash
script="python $HOME/code/wtccc/impute/genotypes.py"
ebi=/broad/compbio/lward/GWAS/ebi_ega_downloads
base=/broad/compbio/aksarkar/projects/impute
diseases=(BD CAD CD HT RA T1D T2D)
chr=({1..22} 24)
for i in ${chr[@]}
do
    for d in ${diseases[@]}
    do
        in=$ebi/WTCCC_7_Diseases/$d/Affymetrix500K/Chiamo/OXSTATS_FORMAT
        out=$base/genotypes/$(echo $d | tr '[A-Z]' '[a-z]')
        echo "zcat $in/${d}_$(printf '%02d' $i)_chiamo.gz | $script $base/affy/chr$i.bed | gzip >$out/$i.gen.gz"
    done
    echo "zcat $ebi/1958BC/WTCCC1_1958BC/Affymetrix500K/Chiamo/Oxstat_format/58C_$(printf '%02d' $i)_chiamo.gz | $script $base/affy/chr$i.bed | gzip >$base/genotypes/1958bc/$i.gen.gz"
    echo "zcat $ebi/NBS/WTCCC1/Affymetrix500K/Chiamo/Oxstat_format/NBS_$(printf '%02d' $i)_chiamo.gz | $script $base/affy/chr$i.bed | gzip >$base/genotypes/nbs/$i.gen.gz"
done
