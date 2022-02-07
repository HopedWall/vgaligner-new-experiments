# If line starts with >, it will be the name for
# the file containing the path
/^>/ {
    OUT=prefix substr($0,2) ".fa"
}

# Otherwise, the line represents the sequence itself,
# so store in in another the file having the name
# defined above
OUT {
    print $0 >OUT
}
