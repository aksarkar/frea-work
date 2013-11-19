# Construct a Plink map file for only well-imputed SNPs
#
# Expects IMPUTE2 info file on stdin, writes map file on stdou.
#
# Author: Abhishek Sarkar <aksarkar@mit.edu>
BEGIN {
    if (chr == "") {
        chr = 0
    }
}

$2 == "." {
    $2 = chr ":" $3
}

!/snp_id/ && $5 > .6 {
    print chr, $2, 0, $3
}
