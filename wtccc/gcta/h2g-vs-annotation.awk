/G0/ {
    f = FILENAME
    sub(/.hsq/, "", f)
    print f, $2, $4
    exit 0
}
