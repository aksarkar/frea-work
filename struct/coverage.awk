{
    s += $3 - $2
}
END {
    print s, s/3e9
}
