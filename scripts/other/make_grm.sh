#!/bin/bash

cd ~/bin/gemma/bin

BIMBAM=$1
PHENO=$2
GRM=$3

./gemma \
    -g $BIMBAM \
    -p $PHENO \
    -gk 2 \
    -o $GRM
