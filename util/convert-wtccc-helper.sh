#!/bin/bash
#BSUB -J convert-wtccc[1-22]
#BSUB -R rusage[argon_io=2]
#BSUB -o convert.log
#BSUB -q compbio-week
awk -vOFS='\t' -vc=$LSB_JOBINDEX 'NR > 1 {print "chr"c, $3 - 1, $3, $2, $16}' /broad/compbio/lward/GWAS/ebi_ega_downloads/WTCCC_7_Diseases/WTCCC_summary_data/7_Diseases/$DISEASE/combined_controls/snptest_${DISEASE}_$(printf "%02d" $LSB_JOBINDEX).txt >/broad/hptmp/aksarkar/wtccc/$DISEASE-$(printf "%02d" $LSB_JOBINDEX).bed
