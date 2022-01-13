#!/usr/bin/env Rscript

library(dplyr)
source("general_fncs.R")

args      = commandArgs(trailingOnly=TRUE)
gex_file  = args[1]
samp_file = args[2]
out_file  = args[3]


################################
############# MAIN #############
################################

## read in RSEM TPM values for genes
gex = read.delim(gex_file, sep="\t", header=T, stringsAsFactors=F)
samps = readLines(samp_file)

gex = gex %>% select(c("gene_id",samps))

dim.gex = dim(gex)

## filter genes that don't have a TPM in at least 5 individuals
gex.clean = gex %>% filter(apply(gex[,samps], 1, function(x) sum(x >= 0.1) > 4))

## quantile normalize
mean_norm = quant_norm_mean(gex.clean[,samps])
norm_norm = quant_norm_norm(mean_norm)

colnames(norm_norm) = samps
norm_norm$gene_id = gex.clean$gene_id
    
write.table(norm_norm[,c("gene_id",samps)],file=out_file,sep="\t",quote=F,row.names=F,col.names=T)

