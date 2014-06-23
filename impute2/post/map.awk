# Construct a Plink map file for only well-imputed SNPs
#
# Expects IMPUTE2 info file on stdin, writes map file on stdout.
#
# Author: Abhishek Sarkar <aksarkar@mit.edu>
NR == 1 {
    split(FILENAME, a, ".")
    chr = a[1]
}

$2 == "." {
    $2 = chr ":" $3
}

NR > 1 && $4 >= .01 && $4 <= .99 && $5 >= .6 {
    print chr, $2, 0, $3
}
