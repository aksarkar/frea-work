#!/bin/bash
tr [A-Z] [a-z] <$HOME/gwas/meta/celltypes.txt | parallel "zcat $HOME/hp/chromhmm/{}.bed.gz | grep Promoter | intersectBed -a $1 -b stdin -c | cut -f6 >{}.annot"
echo "rsid,p,gm12878,h1,hepg2,hmec,hsmm,huvec,k562,nhek,nhlf"
zcat $1 | cut -f4,5 | sed -e "s/\t/,/" | paste -d, - *.annot | sort -t, -g -k2 | grep -v ",,"
rm *.annot
