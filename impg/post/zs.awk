BEGIN {OFS = "\t"}
! /SNP/ && $6 > .6 {print "chr"c, $2-1, $2, $1, $5}
