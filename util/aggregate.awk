BEGIN {
    OFS="\t"
}
{
    if ($4 ~ /Promoter/) {$4 = "promoter"; out = "promoter.bed"}
    else if ($4 ~ /Enhancer/) {out = "enhancer.bed"}
    else if ($4 ~ /Insulator/) {out = "insulator.bed"}
    else if ($4 ~ /Repressed/) {out = "repressed.bed"}
    else if ($4 ~ /Txn/) {out = "transcribed.bed"}
    else {out = "other.bed"}
    print >out
}
