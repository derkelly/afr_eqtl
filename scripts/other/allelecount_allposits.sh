#!/bin/bash

TXTIN=$1
POSITS=$2
TXTOUT=$3

tail -n +2 $TXTIN | \
    awk '{ print $1 "\t" $2-1 "\t" $2 "\t" $4 "\t" $5 "\t" $6 "\t" $7 }' | \
    bedtools intersect -a $POSITS -b stdin -wb -loj | \
    awk '/^[^#]/ { if ($10 == ".") print $3 "\t" $4 "\t" $5 "\t0,0"; else print $3 "\t" $4 "\t" $5 "\t" $15 "," $16; }' | \
    bgzip > $TXTOUT
