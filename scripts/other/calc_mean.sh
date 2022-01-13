#!/bin/bash
#
# This script calculates the mean of the column given, excluding negative values
FILE=$1
COL=$2

tail -n+2 $FILE | cut -f$2 | awk '{ if ($1 == "-nan") ; else if ($1 < 0 ) n++; else sum += $1; n++; } END { if (n > 0) print sum / n; }' 
