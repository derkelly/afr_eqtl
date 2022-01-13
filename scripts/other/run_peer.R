#!/usr/bin/env Rscript

library(peer)
library(dplyr)

source("general_funcs.R")

args = commandArgs(trailingOnly=TRUE)
data_file = args[1]
samp_file = args[2]
thresh    = args[3]
n_thresh  = args[4]
n_factors = args[5]
outfile = args[6]


################################
############# MAIN #############
################################

## read in data
gex = read.delim(data_file,sep="\t",header=T,stringsAsFactors=F)

## read in list of samples to analyze
if (gex %>% select(starts_with("Sample_")) %>% length){
    samps = paste("Sample", readLines(samp_file), sep="_")
} else{
    samps = readLines(samp_file)
}

## check if sample names begin with 'Sample'

## find which variants are above the threshold in at least 'n_thresh' individuals
above_thresh = apply(gex %>% select(samps),1,function(x) sum(x > thresh)) >= n_thresh

## pull out the genes that are above the threshold
gex.clean = gex[above_thresh,samps]

## normalize matrix
gex.norm.mat = gex.clean %>% 
	     quant_norm_mean %>%	
	     quant_norm_norm %>%
	     as.matrix

## run peer
factors = run_peer(t(gex.norm.mat),n_factors)

## place factors into data.frame
factors.df = data.frame(factors)
factors.df$Sample = samps
factors.df = factors.df[,c(ncol(factors.df),1:(ncol(factors.df)-1))]

## write output
write.table(factors.df,outfile,quote=F,sep="\t",row.names=F,col.names=F)
