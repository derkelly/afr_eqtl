#!/bin/bash

BAM=$1
BED=$2

samtools view -h $BAM | \
grep -P "vW|@" | \
samtools view -S -b - | \
bedtools bamtobed -i stdin -tag vW | \
bgzip > $BED
