NR==1 {
    for (i=1; i<=NF; i++) {
        ix[$i] = i
    }
}
NR>1 {
    print $ix[c1], $ix[c2], $ix[c3]
}
