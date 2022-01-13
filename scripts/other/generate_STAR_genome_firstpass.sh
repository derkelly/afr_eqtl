#!/bin/bash
#
# This is generate_STAR_genome_firstpass.sh
#
# This script generates a genome index for use with STAR (v
# 2.5.3a). The genome annotation used will be the GTEx gencode v19
# index, which includes regions masked due to ambiguity from
# overlapping, anti-sense exons, which can confound quantification
# when using unstranded RNA-seq data.

GENDIR="/project/tishkofflab/rna_redo/afr_wbrna/data/genomes"

STAR --runThreadN 8 \
    --runMode genomeGenerate \
    --genomeDir ${GENDIR}/star/firstpass \
    --genomeFastaFiles ${GENDIR}/Homo_sapiens_assembly19.fasta \
    --sjdbGTFfile ${GENDIR}/gencode/gencode.v19.transcripts.patched_contigs.gtf \
    --sjdbOverhang 100
