#!/bin/bash

POP=$1
INDIR=$2
OUTDIR=$3

ZFILES=($( ls ${INDIR}/*${POP}*weir.fst.zscore ))

OUTFILE=${OUTDIR}/${POP}_dstat.txt.gz

paste <( cut -f1,2 ${ZFILES[0]} ) \
      <( paste <( cut -f3 ${ZFILES[0]} ) \
               <( cut -f3 ${ZFILES[1]} ) \
               <( cut -f3 ${ZFILES[2]} ) \
               <( cut -f3 ${ZFILES[3]} ) \
               <( cut -f3 ${ZFILES[4]} ) \
               <( cut -f3 ${ZFILES[5]} ) \
               <( cut -f3 ${ZFILES[6]} ) \
               <( cut -f3 ${ZFILES[7]} ) | \
         awk '{ print $1 + $2 + $3 + $4 + $5 + $6 + $7 + $8 }' ) | \
         bgzip > $OUTFILE
