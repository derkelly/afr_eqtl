#!/bin/bash

module load R-3.2.2

SQTL=$1
FDR=$2
OUT=$3


Rscript get_sig_sqtl.R $SQTL $FDR | bgzip > $OUT
