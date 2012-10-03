#!/bin/bash
#BSUB -J count-maf
#BSUB -q compbio-week
#BSUB -o fix.log
set -e
base=/broad/compbio/lward/GWAS/ebi_ega_downloads
# in=$(find $base -type d -name "imputed" | sed -n "$LSB_JOBINDEX p")
in=/broad/compbio/lward/GWAS/ebi_ega_downloads/WTCCC_7_Diseases/T1D/imputed
out=$(echo $in | sed -e "s#$base##" -e "s#.##" -e "s#/#_#g" | tr [A-Z] [a-z])
for g in $(find $in -name "*.txt.gz" | sort)
do
    chr=$(echo $g | sed -re 's/.*chr([0-9]+).*/\1/')
    zcat $g | grep -v "0ustar" | python $HOME/code/wtccc/maf/maf.py $chr >>$out.tmp
done
~lward/bin/liftOver/liftOver $out.tmp ~lward/bin/liftOver/hg17ToHg19.over.chain $out $out.log
gzip $out
rm $out.tmp
