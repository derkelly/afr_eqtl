#!/bin/bash

HOMEDIR="/home/derkelly"
PROJDIR="/project/tishkofflab/rna_redo/afr_wbrna"
FASTA=${PROJDIR}/data/genomes/Homo_sapiens_assembly19.fasta
TEMPBED=${HOMEDIR}/.tmp/tmp.bed
MERGEBED=${HOMEDIR}/.tmp/tmp.merge.bed
TEMPFA=${HOMEDIR}/.tmp/tmp.fa
OUT=${PROJDIR}/data/misc/secondpass.featureCounts.gc_content.txt

echo -e "gene_id\ttotal_bp\tgc_bp" > $OUT

# read over line, split into chr and bounds for bed file, extract sequence, and calculate gc content.
while read LINE;
do

    GENE=$(echo $LINE | cut -d' ' -f1)
    readarray CHRS < <(echo $LINE | cut -d' ' -f2)
    START=$(echo $LINE | cut -d' ' -f3)
    END=$(echo $LINE | cut -d' ' -f4)
    STRAND=$(echo $LINE | cut -d' ' -f5)
    paste <(echo $CHRS | sed 's/chr//g' | tr ";" "\n") \
	  <(echo $START | tr ";" "\n" | awk '{print $1-1}') \
	  <(echo $END | tr ";" "\n") \
	  <(echo $STRAND | tr ";" "\n") > $TEMPBED

    bedtools merge -i $TEMPBED > $MERGEBED
    bedtools getfasta -s -fi $FASTA -bed $MERGEBED -fo $TEMPFA
    
    SEQ=$(grep "^[^>]" $TEMPFA | tr -d '\n')
    # SEQ_LEN=$(echo "${#SEQ}")

    GC=$(echo $SEQ | sed 's/[AaTt]//g')
    # GC_LEN=$(echo "${#SEQ}")
    
    echo -e "${GENE}\t${#SEQ}\t${#GC}" >> $OUT
done < /dev/stdin
