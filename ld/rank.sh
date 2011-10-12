#!/bin/bash
#BSUB -J ld
#BSUB -o ld.%J
#BSUB -q compbio-week
#BSUB -R rusage[iodine_io=3]

export TMPDIR=/broad/shptmp/aksarkar

# Ranklist of GWAS SNPs
ranklist=ranklist
sort -g -k5 $GWAS >$ranklist

# Independent enhancers
enh=enh

inld=inld
work=work
# while [[ -s $ranklist ]]
# do
    # Take the next SNP on the ranklist
    head -n1 $ranklist | cut -f1-4 >$work
    chr=$(sed -re "s/chr([[:digit:]]+).*/\1/" $work)
    pos=$(cut -f3 $work)
    # Expand top SNP using LD (R^2 > .8)
    zcat $HOME/hp/hapmap/ld_chr${chr}_CEU.txt.gz | awk -v OFS='\t' -v pos=$pos -v chr=$chr -v thresh=.8 --file $HOME/code/ld/ld.awk >$inld
    # Find all such SNPs that fall in enhancers
    bzcat $HOME/hp/chromhmm/hg18/k562-states.bed.bz2 | grep Enhancer | intersectBed -a $inld -b stdin | intersectBed -a stdin -b $ranklist >$work
    # Find representative in enhancer
    head -n1 $work
    # Expand representative
    pos=$(tail -n1 $enh | cut -f3)
    zcat $HOME/hp/hapmap/ld_chr${chr}_CEU.txt.gz | awk OFS='\t' pos={} thresh=.2 --file $HOME/code/ld/ld.awk >$inld
#     # Remove everything in LD with representative
#     intersectBed -v -a $ranklist -b $inld >$work
#     mv $work $ranklist
# done
# rm $ranklist $inld $work
# mv $enh enhancers.out
