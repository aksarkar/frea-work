# -*- makefile -*-
# Makefile for pileup visualization pipeline
# Author: Abhishek Sarkar <aksarkar@mit.edu>

MARKERS=
PFLAGS=1000
RFLAGS=-t .2

all: pileup.pdf vis.png

.PHONY: all

pileup.pdf: pileup expected $(HOME)/code/enr/plot/pileup.R
	R --quiet --vanilla --args pileup expected <$(HOME)/code/enr/plot/pileup.R

vis.png: elements $(HOME)/code/enr/plot/vis.R
	R --quiet --vanilla --args elements <$(HOME)/code/enr/plot/vis.R

chrom: regions.bed exons known_lincRNAs
	cat exons known_lincRNAs | \
    bedtools sort | \
    bedtools merge | \
    bedtools subtract -a $(HOME)/hp/chromhmm/gm12878.bed.gz -b stdin | \
    bedtools intersect -a regions.bed -b stdin -wb | \
    awk -vOFS='\t' '{print $$1, $$2, $$3, $$8, $$4}' | \
    sed -rf $(HOME)/code/util/aggregate.sed >$@

clean:
	rm -f chrom 
	rm -f elements
	rm -f exons 
	rm -f known_lincRNAs
	rm -f pileup
	rm -f pileup.pdf
	rm -f regions.bed
	rm -f vis.png

elements: chrom exons known_lincRNAs
	cat chrom exons known_lincRNAs | \
    python $(HOME)/code/ld/torel.py >$@

exons: regions.bed
	zcat $(HOME)/hp/gencode/gencode.v10.annotation.gtf.gz | \
    awk '$$3 == "exon" && /protein_coding/ && /KNOWN/'  | \
    bedtools intersect -a regions.bed -b stdin -wb | \
    awk -vOFS='\t' '{print $$1, $$2, $$3, "exon", $$4}' >$@

known_lincRNAs: regions.bed
	zcat $(HOME)/hp/gencode/gencode.v10.annotation.gtf.gz | \
    awk -vOFS='\t' '$$3 == "exon" && /lincRNA/ && /KNOWN/' | \
    bedtools intersect -a regions.bed -b stdin -wb | \
    awk -vOFS='\t' '{print $$1, $$2, $$3, "known_lincRNA", $$4}' >$@

pileup: chrom exons known_lincRNAs
	cat chrom exons known_lincRNAs | \
    python $(HOME)/code/ld/torel.py | \
    python $(HOME)/code/util/pileup.py $(BINSIZE) >pileup

expected:
	false
