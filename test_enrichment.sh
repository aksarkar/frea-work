#!/bin/bash
#BSUB -J gwas_chromatin_enrichment[1-4216]
#BSUB -o out/python3.o%J.%I
#BSUB -P compbiofolk
#BSUB -q compbio-week
#BSUB -R "rusage[mem=5, iodine_io=2]"
basepath='/seq/compbio-hp/GWAS/'
python3 "$basepath/enrichment/enrichment.py" $(sed -nre "$LSB_JOBINDEX s#./.*/#$basepath/analyses/#p" "$basepath/analyses.txt")
