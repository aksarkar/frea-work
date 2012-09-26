#!/bin/bash
#BSUB -J postprocess
#BSUB -R rusage[argon_io=3]
#BSUB -o postprocess.log
#BSUB -q compbio-week
cat cases.?? | sort -t' ' -k1 >cases
cat controls_1958bc.?? | sort -t' ' -k1 >controls_1958c
cat controls_nbs.?? | sort -t' ' -k1 >controls_nbs
python $HOME/code/wtccc/merge.py cases controls_1958c controls_nbs >merged
rm cases* controls*
