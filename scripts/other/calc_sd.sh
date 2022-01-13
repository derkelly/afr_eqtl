#!/bin/bash
#
# This script calculates the mean of the column given, excluding negative values
# POP=$1
# COL=$2

MEAN=$1

awk -v mean="$MEAN" '{ if ($1 == "-nan") ; else if ($1 < 0) sum += $mean^2; else sum += ($1-$mean)^2; n++; } END { if (n > 0) print sqrt(sum / n); }'
