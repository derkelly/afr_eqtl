#!/bin/bash
#
# This is map_chunk_STAR_firstpass.v2.7.sh
#
# This script maps reads from multiple samples at a time using
# STAR. This takes a file with samples to run, a directory for where
# to find reads, and an output directory. I do not believe any special
# parameters are required, aside from the STAR genome generated with
# the masked GTEx gene reference.

SAMPFILE=$1
OUTDIR=$2

DATADIR="/project/tishkofflab/rna_redo/afr_wbrna/data"

READDIR=$DATADIR/samples
STARGENO=$DATADIR/genomes/star/firstpass_wasp
VCF=$DATADIR/geno/5M.imputed.dose.biallelic.dummy.vcf

readarray -t SAMPLES<$SAMPFILE

# map all samples in a given chunk
for SAMPLE in ${SAMPLES[@]};
do

    SAMPOUT=${OUTDIR}/${SAMPLE}
    PREF=${SAMPOUT}/firstpass_wasp
    
    mkdir -p $SAMPOUT

    FORWARD=${READDIR}/${SAMPLE}/forward.noribo.fa.gz
    REVERSE=${READDIR}/${SAMPLE}/reverse.noribo.fa.gz

    # echo $PREF
    # echo $STARGENO
    # echo $FORWARD
    # echo $REVERSE
    # echo $VCF
    
    ~/bin/STAR-2.7.0f/source/STAR \
    	--outFileNamePrefix $PREF \
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
    	--quantMode TranscriptomeSAM \
    	--varVCFfile $VCF \
    	--waspOutputMode SAMtag \
    	--outSAMtype BAM SortedByCoordinate \
    	--limitBAMsortRAM 10000000000

done

~/bin/STAR-2.7.0f/source/STAR --genomeDir $STARGENO --genomeLoad Remove
