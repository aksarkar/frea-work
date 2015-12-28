BEGIN {
    while ((getline < (PREFIX "/prevalences.txt")) > 0) {
        phenos[$1] = $2
    }
    for (p in phenos) {
        for (i = 1; i < 23; i++) {
            print "bash ~/code/wtccc/gcta/hsq_vs_n.sh", PREFIX "/" p, phenos[p], "observed", i >"joblist."p
            print "bash ~/code/wtccc/gcta/hsq_vs_n.sh", PREFIX "/" p, phenos[p], "random", i >"joblist."p
        }
    }
}
