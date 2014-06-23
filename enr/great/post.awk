BEGIN {
    FS="\t"
    OFS="\t"
}

NR == 1 {
    f=FILENAME
    sub(/^.*\//, "", f)
    sub(/\..*$/, "", f)
    print "cluster", "term", "fold", "p", "q", "genes"
}

NR == 5 {
    sub(/^# /, "", $0)
    for (i = 1; i < NF; i++) {
        header[$i] = i
    }
}

(NR > 5 && 
 $header["Ontology"] == "GO Biological Process" && 
 $header["HyperFdrQ"] < 1e-4 && 
 $header["RegionFoldEnrich"] > 1) {
    print f, $header["Desc"], $header["RegionFoldEnrich"], $header["HyperP"], $header["HyperFdrQ"], $header["FgGeneNames"]
}
