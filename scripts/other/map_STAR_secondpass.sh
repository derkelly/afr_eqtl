#!/bin/bash
#
# This is map_STAR_second.sh
#
# This script maps reads from a single sample at a time.

SAMPDIR=$1

STARGENO="/project/tishkofflab/rna_redo/afr_wbrna/data/genomes/star/secondpass"

SAMPOUT=$SAMPDIR
PREF=${SAMPDIR}/secondpass

mkdir -p $SAMPDIR

FORWARD=${SAMPDIR}/forward.noribo.fa.gz
REVERSE=${SAMPDIR}/reverse.noribo.fa.gz

STAR --outFileNamePrefix $PREF \
    --genomeLoad LoadAndKeep \
    --genomeDir $STARGENO \
    --outSAMunmapped Within \
    --readFilesCommand zcat \
    --readFilesIn $FORWARD $REVERSE \
    --runThreadN 16 \
    --outFilterType BySJout \
    --outSAMattributes NH HI AS NM MD \
    --outFilterMultimapNmax 20 \
    --outFilterMismatchNmax 999 \
    --outFilterMismatchNoverReadLmax 0.04 \
    --alignIntronMin 20 \
    --alignIntronMax 1000000 \
    --alignMatesGapMax 1000000 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 1 \
    --sjdbScore 1 \
    --quantMode TranscriptomeSAM

STAR --genomeDir $STARGENO --genomeLoad Remove
