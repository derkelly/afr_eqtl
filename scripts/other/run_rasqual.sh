#!/bin/bash

FEATURES=$1
VCF=$2
Y=$3
K=$4
X=$5
N=$6
WIN=$7
SCRIPT=$8
OUTDIR=$9

MYTMPDIR=$(mktemp -d)
trap "rm -rf $MYTMPDIR" EXIT

#j=1 # index variable
bsub j=$LSB_JOBINDEX

