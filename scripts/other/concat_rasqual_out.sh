#!/bin/bash
#
# This is concat_rasqual_out.sh
#
# This script takes a directory and produces a 'bgzip'-ed and
# 'tabix'-ified concatenation of all RASQUAL results in that directory

INDIR=$1
OUTFILE=$2

cat $INDIR/* | \
    sort -k1,1 -k4,4n | \
    bgzip > $OUTFILE

tabix -s 1 -b 4 -e 4 $OUTFILE
