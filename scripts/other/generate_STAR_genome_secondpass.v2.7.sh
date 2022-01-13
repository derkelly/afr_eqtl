#!/bin/bash
#
# This is generate_STAR_genome_secondpass.v2.7sh
#
# This script generates a genome index for use with STAR (v
# 2.7.0f). The genome annotation used will be the GTEx gencode v19
# index, which includes regions masked due to ambiguity from
# overlapping, anti-sense exons, which can confound quantification
# when using unstranded RNA-seq data. In addition, filtered splice
# junctions will be used for the second pass of mapping.

DATADIR="/project/tishkofflab/rna_redo/afr_wbrna/data"
GENDIR=${DATADIR}/genomes
OUTDIR=$GENDIR/star/secondpass_wasp

STAR --runThreadN 8 \
    --runMode genomeGenerate \
    --genomeDir $OUTDIR \
    --genomeFastaFiles ${GENDIR}/Homo_sapiens_assembly19.fasta \
    --sjdbGTFfile ${GENDIR}/gencode/gencode.v19.transcripts.patched_contigs.gtf \
    --sjdbOverhang 100 \
    --sjdbFileChrStartEnd ${DATADIR}/samples/Sample_ALL/firstpassSJ.all.canon.1reads.2samp.out.tab \
    --limitSjdbInsertNsj 4000000
