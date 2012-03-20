#!/bin/bash
#BSUB -J control[1-500]
#BSUB -o log
#BSUB -q compbio-week
#BSUB -R rusage[mem=2,iodine_io=3]
t="$LSB_JOBID.$LSB_JOBINDEX"
mkdir $t
cd $t
ln -s $HOME/code/util/pileup.make Makefile
touch acc
for i in {0..99}
do
    zcat $HOME/hp/t1d/data/hg19/t1d.bed.gz | \
        python3 $HOME/code/util/sample.py 30000 | \
        gzip >markers.bed.gz
    make -s pileup MARKERS=markers.bed.gz RFLAGS="-s 200000"
    python $HOME/code/util/acc.py acc pileup >temp
    mv temp acc
    rm regions.bed
done
