#!/bin/bash

VCF=$1
OUT=$2

plink --vcf $VCF \
      --recode bimbam \
      --out $OUT
