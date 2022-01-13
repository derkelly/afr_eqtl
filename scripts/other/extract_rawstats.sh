#!/bin/bash
#
# This is extract_rawstats.sh
#
# This script extracts the feature ID, SNP ID, ref allele, alt allele,
# chi-square statistic, and direction of the effect size from RASQUAL
# and outputs the ID, reftput for further processing.

INDIR=$1
N=$2
OUT=$3

cat ${INDIR}/* | \
    awk '{ if ($12 > 0.5) print $1 ":" $2 "\t" $5 "\t" $6 "\t" $11 "\t+\t" $N;
    	   else print $1 ":" $2 "\t" $5 "\t" $6 "\t" $11 "\t-\t" $N}' | \
	       Rscript -e 'data=read.delim("stdin",sep="\t",header=F,stringsAsFactors=F);\
	       data$chisq_p=pchisq(data[,4],df=1,lower.tail=F); \
	       write.table(data[,c(1:3,7,5,6)],sep="\t",col.names=F,row.names=F,quote=F)'> $OUT
