#!/bin/bash

TESTREGION=$1
VCF=$2
Y=$3
K=$4
X=$5
N=$6
j=$7
NTESTVARS=$8
NFEATVARS=$9
MERGESTART=${10}
MERGEEND=${11}
ID=${12}
OUT=${13}

tabix $VCF $TESTREGION | ~/bin/rasqual/bin/rasqual -y $Y \
    -k $K \
    -x $X \
    -n $N \
    -j $j \
    -l $NTESTVARS \
    -m $NFEATVARS \
    -s $MERGESTART \
    -e $MERGEEND \
    -r \
    -f $ID > $OUT
