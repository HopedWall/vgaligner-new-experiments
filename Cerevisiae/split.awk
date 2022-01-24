/^>/ {
    OUT=substr($0,2) ".fa"
}

OUT {
    print >OUT
}