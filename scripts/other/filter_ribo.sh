#!/bin/bash
#
# This is filter_ribo.sh
#
# This function takes a file with read ids to exclude and removes them
# from a fasta file. This assumes the ids to remove are in the same
# relative order in each file.

FASTA=$1
RIBO=$2
OUT=$3

faSomeRecords -exclude \
$FASTA $RIBO $OUT

bgzip $OUT
