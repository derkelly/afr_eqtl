#!/bin/bash
#
# This is filter_starout.sh
#
# This function filters all output from the second pass of STAR

## input includes the main directory for output and the name of the sample
MAINDIR=$1
SAMPLE=$2

# output directory
SAMPDIR=${MAINDIR}/${SAMPLE}

## set input file names
GENOIN=${SAMPDIR}/secondpassAligned.out.sam
TRXOIN=${SAMPDIR}/secondpassAligned.toTranscriptome.out.bam

## set temporary file names
GENOFLT=${SAMPDIR}/secondpassAligned.out.255.bam
GENOSRT=${SAMPDIR}/secondpassAligned.out.255.sort.bam

## set output file names
GENOOUT=${SAMPDIR}/secondpassAligned.out.255.sort.rg.bam
TRXOOUT=${SAMPDIR}/secondpassAligned.toTranscriptome.out.255.bam

## filter files
samtools view -bhq 255 $GENOIN \
    -o $GENOFLT

samtools view -bhq 255 $TRXOIN \
    -o $TRXOOUT

## sort file
samtools sort -@ 8 -m 1G -o $GENOSRT $GENOFLT

## add read group
module load java-sdk-1.8.0
java -jar ~/bin/picard.jar AddOrReplaceReadGroups \
    I=${GENOSRT} \
    O=${GENOOUT} \
    RGLB=LaneX \
    RGPL=illumina \
    RGPU=NONE \
    RGSM=${SAMPLE}

## remove unnecessary temporary files
rm $GENOFLT
rm $GENOSRT
