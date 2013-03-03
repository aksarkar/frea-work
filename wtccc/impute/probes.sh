# Lift over Affy probe positions
#!/bin/bash
cd /broad/compbio/lward/GWAS/ebi_ega_downloads/WTCCC_7_Diseases/T1D/Affymetrix500K/Chiamo/OXSTATS_FORMAT
for i in (1..24)
do 
    f=T1D_$(printf "%02d" $i)_chiamo.gz
    if [[ -f $f ]] 
    then
        zcat $f | awk -vOFS='\t' -vc=$i '{print "chr"c, $3-1, $3, $1}' >>probes-hg17.bed
    fi
done
~lward/bin/liftOver/liftOver probes.bed ~lward/bin/liftOver/hg17ToHg19.over.chain out log
bedtools sort <out >probes.bed
