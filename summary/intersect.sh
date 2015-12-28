#!/bin/bash
ln -s /broad/compbio/aksarkar/projects/gwas/wtccc1/EC21/results/impg/out/*.bed.gz .
parallel --joblog 0.log -n2 "bedtools intersect -wa -sorted -nobuf -a {1} -b {2} | bedtools sort >{#}.0" ::: *.bed.gz
parallel --joblog 1.log -n2 "bedtools intersect -wa -sorted -nobuf -a {1} -b {2} | bedtools sort >{#}.1" ::: *.0
parallel --joblog 2.log -n2 "bedtools intersect -wa -sorted -nobuf -a {1} -b {2} | bedtools sort >{#}.2" ::: *.1
gzip <1.2 >all-studies.bed.gz
