#!/bin/bash
#BSUB -J convert-wtccc[1-1]
#BSUB -o /dev/null
#BSUB -q compbio-week
diseases=(BD CAD CD HT RA T1D T2D)
DISEASE=${diseases[$LSB_JOBINDEX - 1]} bsub <$HOME/code/util/convert-wtccc-helper.sh
