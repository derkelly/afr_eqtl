#!/bin/bash

VCF=$1
SAMPS=$2
OUT=$3

vcftools --gzvcf $VCF \
	 --keep $SAMPS \
	 --freq \
	 --out $OUT
