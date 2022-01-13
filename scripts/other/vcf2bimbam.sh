#!/bin/bash

VCF=$1
OUT=$2

zcat $VCF | \
    cut -f3,4,5,10-171 | \
    sed 's/.\|.\://g
awk '{for(i=1;i<=NF;i++)x+=$i;print x}'

paste <(zcat $VCF | \
	    grep -vE "^#") \
      <(zcat $VCF | \
	    grep -vE "^#" | \
	    perl -nle 'print $1." " ile /\d+\|\d+:([0-9\.]+)/g' )
