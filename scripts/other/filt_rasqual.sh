#!/bin/bash

IN=$1
OUT=$2

zcat $IN | \
    awk '{ if ($13 < 0.1 && $14 > 0.25 && $14 < 0.75 && $24 > 0.9 && $25 > 0.9) print $0 }' | \
    bgzip > $OUT
