#!/bin/bash
#
# This is run_gemma.sh
#
# This script runs gemma for a given phenotype file using genotype
# data in the current directory. This is designed to run eQTL mapping

# This function calls GEMMA, sending the output to the file OUTTMP to
# be processed
GEMMA_CALL(){

    GENOTMP=$1
    PHENOTMP=$2
    GRM=$3
    OUTTMP=$4

    ~/bin/gemma/bin/gemma \
    	-g $GENOTMP \
    	-p $PHENOTMP \
    	-k $GRM \
    	-lmm 4 \
    	-o $OUTTMP
}


# This function generates temporary files necessary for running GEMMA,
# namely the GENOTMP genotype file and the PHENOTMP phenotype file
MAKE_FILES(){

    LINE=$1
    GENO=$2
    GENOTMP=$3
    PHENOTMP=$4

    # write current phenotype to file
    echo $LINE | \
	cut -d' ' -f1-5 --complement | \
	sed 's/ /\n/g' > $PHENOTMP

    # get the range within 100kb of the TSS
    SNPRANGE=$(echo "$LINE" | \
	awk '{ if ($4=="+") print $1 "\t" $2-1 "\t" $2; else print $1 "\t" $3-1 "\t" $3 }' | \
	bedtools slop -i stdin -g $CHRSIZE -b 100000 | \
	awk '{ print $1 ":" $2 "-" $3 }')

    # extract relevant genotypes
    tabix $GENO $SNPRANGE | \
	cut -f1,2 --complement > $GENOTMP
}

##################################################
###################### MAIN ######################
##################################################

# take as input the phenotype file to run on
INPATH=$1

# get the filename and the number of peer factors
INFILE=$(basename $INPATH)
# PEER=$( echo $INPATH | perl -nE 'print $1."\n" if /.*?(peer\d+)/' )

# make a termporary directory
TMPDIR="tmp"
mkdir -p $TMPDIR

# set files and paths that will be necessary
GENO="5M.imputed.dose.biallelic.doses.txt.gz"
GRM="5M.imputed.dose.biallelic.doses.grm.sXX.txt"
OUTTMP=$INFILE
GENOTMP=$TMPDIR/$INFILE.geno.tmp.txt
PHENOTMP=$TMPDIR/$INFILE.pheno.tmp.txt
CHRSIZE="/local1/derek/data/genomes/hg19.nochr.chrom.sizes"

# OUTDIR=results/$PEER
# mkdir -p $OUTDIR

OUTDIR=conditional_map
mkdir -p $OUTDIR

# for each line:
# 1. extract the phenotype information
# 2. make the name for the output file
# 3. make all necessary temporary files, including a temporary genotype and phenotype file
# 4. call gemma
# 5. concatenate the results
while read LINE; do

    if [[ "$LINE" =~ ^[0-9]+ ]]
    then

	PHENO=$(echo "$LINE" | cut -f5)
	GEMMAOUT=$OUTDIR/$OUTTMP."gemma.100kb.lmm4"
	
	MAKE_FILES "$LINE" $GENO $GENOTMP $PHENOTMP
	
	GEMMA_CALL $GENOTMP $PHENOTMP $GRM $OUTTMP
	
	tail -n +2 output/$OUTTMP.assoc.txt | \
    	    sed "s/^/${PHENO}\t/g" >> $GEMMAOUT
	
	rm output/$OUTTMP.assoc.txt
    fi

done <$INPATH
