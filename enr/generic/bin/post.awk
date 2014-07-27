BEGIN {
    FS=","
}

{
    print >$2"."$4".in"
}
