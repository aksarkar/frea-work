BEGIN {
    FS="," 
    OFS=","
}
{
    print >>$2".in"
}
