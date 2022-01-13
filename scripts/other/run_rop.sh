#!/bin/bash

IN=$1
OUTDIR=$2

~/bin/rop/rop.sh -bf -s immune,microbiome $IN $OUTDIR
