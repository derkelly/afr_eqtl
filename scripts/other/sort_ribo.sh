#!/bin/bash
#
# This is sort_ribo.sh
#
# This function takes a file with read ids and sorts them, saving them
# to the 'OUT' file.

RIBO=$1
OUT=$2

sed 's/\./\t/g' $RIBO | \
    sort -k2,2n | \
    sed 's/\t/\./g' | bgzip > $OUT
