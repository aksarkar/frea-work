{
    if (/CDS/ || /3UTR/ || /5UTR/ || /ncRNA/) {
        print >"gencode/"$4".bed"
    }
    else if (/Enh/ || /Prom/ || /Ins/ || /Trans/) {
        split($4, a, /-/)
        if (length(a) == 1) {
            print >$4".bed"
        }
        else {
            print >a[2]"/"a[1]".bed"
        }
    }
    else {
        print >"tfbs/"$4".bed"
    }
}
