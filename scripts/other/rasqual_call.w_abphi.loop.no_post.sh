#!/bin/bash

# This is adapted from rasqual_call.w_abphi.sh. Rather than using the
# 'j' value from the LSB_JOBINDEX it is supplied by the user. I will
# be splitting the feature file into 5 files to run in parallel

FEATURES=$1
VCF=$2
Y=$3
K=$4
X=$5
N=$6
ABPHI=$7
WIN=$8
MYTMPDIR=$9
OUTDIR=${10}
j=${11}

LINE=$(awk "NR==$j" $FEATURES)


ID=$(echo "$LINE" | cut -f1)
STRAND=${LINE: -1}

TMPBED=${MYTMPDIR}/${ID}.bed

# bed-formatted positions of merged exons
MERGEBED=$(paste <(echo "$LINE" | cut -f2 | sed 's/;/\n/g') \
		 <(echo "$LINE" | cut -f3 | sed 's/;/\n/g') \
		 <(echo "$LINE" | cut -f4 | sed 's/;/\n/g') | \
	       bedtools merge)
echo "$MERGEBED" > $TMPBED

# list of feature start positions, separated by commas
CHR=$(echo "$MERGEBED" | cut -f1 | head -n 1)
MERGESTART=$(echo "$MERGEBED" | cut -f2)
MERGEEND=$(echo "$MERGEBED" | cut -f3)

# get the TSS for the test window
TSS=$(echo "$MERGESTART" | head -n 1)

if [ "$STRAND" == "-" ]
then
    TSS=$(echo "$MERGEEND" | tail -n 1)
fi

# change newline characters to commas for RASQUAL
MERGESTART=$(echo $MERGESTART | sed 's/ /,/g')
MERGEEND=$(echo $MERGEEND | sed 's/ /,/g')

# get start and end of test window
TESTSTART=$(($TSS - $WIN))
if [ $TESTSTART -lt 0 ]; 
then 
    TESTSTART=0 
fi
TESTEND=$(($TSS + $WIN))

# get test variants
TESTREGION=${CHR}:${TESTSTART}-${TESTEND}
TESTVARS=$(tabix $VCF $TESTREGION)

if [ ${#TESTVARS} -ne "0" ];
then
    NTESTVARS=$(echo "$TESTVARS" | wc -l)
    NFEATVARS=$(echo "$TESTVARS" | \
		    awk '{ print $1 "\t" $2-1 "\t" $2 }' |
		    bedtools intersect -a stdin -b $TMPBED -u | wc -l)
    OUT=${OUTDIR}/${ID}.rasqual.out

    tabix $VCF $TESTREGION | \
    	~/bin/rasqual/bin/rasqual -ABPHI $ABPHI \
				  --genotype-dosage \
    				  -y $Y \
    				  -k $K \
    				  -x $X \
    				  -n $N \
    				  -j $j \
    				  -l $NTESTVARS \
    				  -m $NFEATVARS \
    				  -s $MERGESTART \
    				  -e $MERGEEND \
				  --n-threads 8 \
				  --no-posterior-update \
    				  -f $ID > $OUT
fi
