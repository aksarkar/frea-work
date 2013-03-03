#!/bin/bash
#BSUB -J prephase[1-22]
#BSUB -q compbio-week
#BSUB -o log
shapeit -G /broad/compbio/aksarkar/wtccc/hg19+dbsnp137/${LSB_JOBINDEX}.gen.txt.gz \
    -M /broad/compbio/aksarkar/genetic_maps_b37/genetic_map_chr${LSB_JOBINDEX}_combined_b37.txt \
    -O $LSB_JOBINDEX.hap.txt $LSB_JOBINDEX.sample
