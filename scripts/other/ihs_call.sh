#!/bin/bash

VCF=$1
MAP=$2
OUT=$3

~/bin/selscan/selscan --ihs --threads 8 --vcf $VCF --map $MAP --out $OUT
