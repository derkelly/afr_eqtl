#!/bin/bash
#
# This script calculates the mean of the column given, excluding negative values

PIPEDIR="/project/tishkofflab/rna_redo/afr_wbrna/pipeline"

COL=$1
MEAN=$2
SD=$3

awk -v col="$COL" -v mean="$MEAN" -v sd="$SD" '{ if ($col == "-nan") print "NA" ; else if ($col < 0) print -$mean/$sd; else sum += ($col-$mean)/$sd; }'
