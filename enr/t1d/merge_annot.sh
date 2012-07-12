#!/bin/bash
echo "rsid,p,gm12878,h1,hepg2,hmec,hsmm,huvec,k562,nhek,nhlf"
zcat $1 | cut -f4,5 | tr '\t' , | paste -d, - *.annot | sort -t, -g -k2 | grep -v ",,"
rm *.annot
