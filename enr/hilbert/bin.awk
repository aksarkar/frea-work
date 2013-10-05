BEGIN {
    OFS="\t"
    binsize = 262144
    chrom = ""
    bin = 0
    count = 0
}

$1 != chrom {
    chrom = $1
    bin = 0
    count = 0
}

$2 > bin + binsize {
    if (count > 0) {
        mid = bin + binsize / 2
        print chrom, mid, mid + 1, ".", count
    }
    while ($2 > bin + binsize) {
        bin += binsize
    }
    count = 0
}

{
    count += 1
}
