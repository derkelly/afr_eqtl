#!/bin/bash


CHR=$1
VCF=$2
HAP=$3
IDX=$4
TAB=$5

   ~/bin/WASP/snp2h5/snp2h5 --chrom $CHR \
         --format vcf \
         --haplotype $HAP \
         --snp_index $IDX \
         --snp_tab $TAB\
         $VCF
