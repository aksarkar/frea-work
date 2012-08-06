#!/bin/bash
#BSUB -J convert-wtccc
#BSUB -R rusage[argon_io=1]
#BSUB -o convert-wtccc.log
#BSUB -q compbio-week
diseases=(BD CAD CD HT RA T1D T2D)
for d in ${diseases[@]}
do
    out=$(echo $d | tr 'A-Z' 'a-z').bed
    for c in {1..22}
    do
        awk -vOFS='\t' -vc=$c 'NR > 1 {print "chr"c, $3 - 1, $3, $2, $16}' /broad/compbio/lward/GWAS/ebi_ega_downloads/WTCCC_7_Diseases/WTCCC_summary_data/7_Diseases/${d}/combined_controls/snptest_${d}_$(printf "%02d" $c).txt >>$d.tmp
    done
    ~lward/bin/liftOver/liftOver $d.tmp ~lward/bin/liftOver/hg17ToHg19.over.chain $out $d.log
    gzip $out
    rm $d.tmp
done
