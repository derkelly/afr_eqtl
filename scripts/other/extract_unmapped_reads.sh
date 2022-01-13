#!/bin/bash

IN=$1
OUT=$2

samtools view -bf 4 $IN | \
    samtools bam2fq - > $OUT
