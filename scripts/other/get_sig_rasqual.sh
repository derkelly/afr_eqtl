#!bin/bash

module load R-3.2.2

INFILE=$1
FDR=$2
OUTFILE=$3

Rscript get_sig_rasqual.R \
	$INFILE \
	$FDR | \
    egrep "^EN" | \
    sort -k1,1 -k4,4n | \
	bgzip -f > $OUTFILE

tabix -s 1 -b 4 -e 4 $OUTFILE
