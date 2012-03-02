# -*- makefile -*-
# Makefile for pileup visualization pipeline
# Author: Abhishek Sarkar <aksarkar@mit.edu>

REGIONS=regions.bed
BINSIZE=1000

out.pdf: pileup
	R --quiet --vanilla --args pileup <$(HOME)/code/enr/plot/pileup.R

chrom: exons known_lincRNAs
	cat exons known_lincRNAs | \
    bedtools sort | \
    bedtools merge | \
    bedtools subtract -a $(HOME)/hp/chromhmm/gm12878.bed.gz -b stdin | \
    bedtools intersect -a $(REGIONS) -b stdin -wb | \
    awk -vOFS='\t' '{print $$4, $$6, $$7, $$8, $$4}' | \
    sed -rf $(HOME)/code/util/aggregate.sed >chrom

exons:
	zcat $(HOME)/hp/gencode/gencode.v10.annotation.gtf.gz | \
    awk '$$3 == "exon" && /protein_coding/ && /KNOWN/'  | \
    bedtools intersect -a $(REGIONS) -b stdin -wb | \
    awk -vOFS='\t' '{print $$5, $$8, $$9, "exon", $$4}' >exons

known_lincRNAs:
	zcat $(HOME)/hp/gencode/gencode.v10.annotation.gtf.gz | \
    awk -vOFS='\t' '$$3 == "exon" && /lincRNA/ && /KNOWN/' | \
    bedtools intersect -a $(REGIONS) -b stdin -wb | \
    awk -vOFS='\t' '{print $$5, $$8, $$9, "known_lincRNA", $$4}' >known_lincRNAs

pileup: chrom exons known_lincRNAs
	cat chrom exons known_lincRNAs | \
    python $(HOME)/code/ld/torel.py | \
    python $(HOME)/code/util/pileup.py $(BINSIZE) >pileup
