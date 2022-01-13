#!/bin/bash

VCF=$1
OUT=$2

zcat $VCF | \
    cut -f1-9 | \
    bgzip > $OUT