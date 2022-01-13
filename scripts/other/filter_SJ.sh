#!/bin/bash
#
# This is filter_SJ.sh
#
# This script concatenates all SJ.out.tab files from the first round
# of STAR mapping and extracts all canonical junctions and all
# junctions with at least 1 unique reads in at least 2 samples.

cat <( cat /project/tishkofflab/rna_redo/afr_wbrna/data/samples/Sample_*/firstpassSJ.out.tab | \
    awk '{ if ($6 == 1) print $0 }' | 
    cut -f1-3 | \
    sort -k1,1 -k2,2n | \
    uniq ) \
    <( cat /project/tishkofflab/rna_redo/afr_wbrna/data/samples/Sample_*/firstpassSJ.out.tab | \
    awk '{ if ($7 > 0) print $0 }' | \
    cut -f1-3 | \
    sort -k1,1 -k2,2n -k3,3n | \
    uniq -c | \
    awk '{ if ($1 > 1) print $0 }' | \
    tr -s ' ' | \
    cut -d' ' -f3 ) | \
	sort -k1,1 -k2,2n -k3,3n | \
	uniq > /project/tishkofflab/rna_redo/afr_wbrna/data/samples/Sample_ALL/firstpassSJ.all.canon.1reads.2samp.out.tab
