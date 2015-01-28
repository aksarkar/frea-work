awk '/narrowPeak/ && !/revoked/ {sub(/cell=/, "", $8); sub(/;/, "", $8); sub(/treatment=/, "", $9); sub (/;/, "", $9); if ($9 != "None") {c=$8"-"$9} else {c=$8} print $1, c}' files.txt | \
    parallel --dry-run -l2 -C' ' 'bedtools intersect -a {1} -b {3} | bedtools sort | bedtools merge -i stdin | gzip >/broad/compbio/aksarkar/data/roadmap/core-features/DHS/ENCODE-{4}.bed.gz'
