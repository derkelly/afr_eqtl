#!/bin/bash
#
# This is sam2bam.sh

IN=$1
OUT=$2

samtools view -Sbh $IN > $OUT
