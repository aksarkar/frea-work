#!/bin/bash
#BSUB -J classifier
#BSUB -o classifier.log
#BSUB -q compbio-week
set -e
python $HOME/code/wtccc/classifier.py merged
python3 $HOME/code/wtccc/avg_roc.py $NSAMPLES $LABEL roc.* >avg.csv
