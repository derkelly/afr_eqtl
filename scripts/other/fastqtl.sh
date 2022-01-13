#!/bin/bash

VCF=$1
BED=$2
COVAR=$3
CHR=$4
OUT=$5

~/bin/FastQTL/bin/fastQTL.static \
    --vcf $VCF \
    --bed $BED \
    --cov $COVAR \
    --window 100000 \
    --permute 1000 10000 \
    --seed 4690 \
    --out $OUT \
    --region $CHR
				 
