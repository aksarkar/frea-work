#!/bin/bash
#BSUB -J classifier
#BSUB -o classifier.log
#BSUB -q compbio-week
set -e
python $HOME/code/wtccc/classifier.py merged
