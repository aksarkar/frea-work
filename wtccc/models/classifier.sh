#!/bin/bash
#BSUB -J classifier
#BSUB -o classifier.log
#BSUB -q compbio-week
set -e
test=/broad/compbio/aksarkar/t1d/models/test.csv
train=/broad/compbio/aksarkar/t1d/models/train.csv
python $HOME/code/wtccc/classifier.py merged $train $test
