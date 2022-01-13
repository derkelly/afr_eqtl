#!/bin/bash
#
# This is extract_rawstats.sh
#
# This script extracts the feature ID, SNP ID, ref allele, alt allele,
# chi-square statistic, and direction of the effect size from RASQUAL
# and outputs the ID, ref allele, alt allele, p-value, and direction.

IN=$1
N=$2
OUT=$3

zcat $IN | \
    awk '{ if ($12 > 0.5) print $1 ":" $2 "\t" $5 "\t" $6 "\t" $11 "\t+\t" '''$N''';
    	   else print $1 ":" $2 "\t" $5 "\t" $6 "\t" $11 "\t-\t" '''$N''' }' | \
	       Rscript -e 'data=read.delim("stdin",sep="\t",header=F,stringsAsFactors=F); data$chisq_p=pchisq(data[,4],df=1,lower.tail=F); data.out=data[,c(1:3,7,5,6)]; colnames(data.out)=c("id","ref","alt","p_val","effect_dir","n"); write.table(data.out,sep="\t",col.names=T,row.names=F,quote=F)' > $OUT
