#!/bin/bash
#BSUB -J normalize
#BSUB -R "rusage[argon_io=3]"
#BSUB -o log
#BSUB -q compbio-week
cat /broad/compbio/aksarkar/ld/1kg/CEU.txt | \
    python $HOME/code/ld/1kg/normalize.py | \
    awk 'NF == 4' >normalized.txt
