#!/bin/bash
#
# This is rsem_call.sh
#
# This script takes a bam file, an RSEM genome, and an output file and
# calls rsem-calculate-expression

BAM=$1
GENO=$2
OUT=$3

rsem-calculate-expression -p 8 \
    --paired-end \
    --alignments \
    --no-bam-output \
    --no-qualities \
    --seed 4690 \
    --paired-end \
    --bam \
    $BAM $GENO $OUT
