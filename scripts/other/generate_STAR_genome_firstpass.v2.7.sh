#!/bin/bash
#
# This is generate_STAR_genome_firstpass.v2.7.sh
#
# This script generates a genome index for use with STAR (v
# 2.7.0f). The genome annotation used will be the GTEx gencode v19
# index, which includes regions masked due to ambiguity from
# overlapping, anti-sense exons, which can confound quantification
# when using unstranded RNA-seq data.

GENDIR="/project/tishkofflab/rna_redo/afr_wbrna/data/genomes"
OUTDIR=$GENDIR/star/firstpass_wasp

mkdir -p $OUTDIR

~/bin/STAR-2.7.0f/source/STAR --runThreadN 8 \
    --runMode genomeGenerate \
    --genomeDir $OUTDIR \
    --genomeFastaFiles ${GENDIR}/Homo_sapiens_assembly19.fasta \
    --sjdbGTFfile ${GENDIR}/gencode/gencode.v19.transcripts.patched_contigs.gtf \
    --sjdbOverhang 100
