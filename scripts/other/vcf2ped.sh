#!/bin/bash

VCF=$1
OUT=$2

zcat $VCF | \
    grep -Pv 'A\tT|T\tA|C\tG|G\tC' | \
    vcftools --vcf - \
    --plink \
    --out $OUT
