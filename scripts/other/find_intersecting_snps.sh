#!/bin/bash

# specify python version
source $HOME/my_python-3.4.2/bin/activate

TAB=$1
IDX=$2
HAP=$3
SAMP=$4
BAM=$5
OUTDIR=$6

python ~/bin/WASP/mapping/find_intersecting_snps.py \
       --is_paired_end \
       --is_sorted \
       --output_dir $OUTDIR \
       --snp_tab $TAB \
       --snp_index $IDX \
       --haplotype $HAP \
       --samples $SAMP \
       $BAM
