#!/bin/bash
# Annotate a list of SNPs with the closest Entrez gene ID for testing
# enrichment of MSigDB pathways
#
# Author: Abhishek Sarkar <aksarkar@mit.edu>
set -e
# map=$(mktemp -p /broad/hptmp/aksarkar)
# zcat /broad/compbio/lward/incoming/entrez/Homo_sapiens.gene_info.gz | \
#     awk -vOFS='\t' '/ENSG[0-9]+/ {match($0, /ENSG[0-9]+/); print substr($0, RSTART, RLENGTH), $2}' | \
#     sort -k1 >$map
map=/broad/hptmp/aksarkar/map
zcat /broad/compbio/lward/incoming/gencode/v13/gencode.v13.annotation.gtf.gz | \
    awk '$3 == "exon" && $14 ~ "\"protein_coding\""' | \
    bedtools closest -a $1 -b stdin -t first | \
    awk -vOFS='\t' '{match($0, /ENSG[0-9]+/); print $1, $2, $3, substr($0, RSTART, RLENGTH), $5}' | \
    sort -k4 | \
    join - $map -14 -t$'\t' -o '1.1 1.2 1.3 2.2 1.5' | \
    bedtools sort
