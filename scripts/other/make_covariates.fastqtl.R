#!user/bin/env Rscript

library(dplyr)
library(tidyr)



args = commandArgs(trailingOnly=TRUE)


## set some useful directories
datadir = "/project/tishkofflab/rna_redo/afr_wbrna/data"
sampdir = file.path(datadir, "samples")
rasqdir = file.path(datadir, "rasqual")
fastdir = file.path(datadir, "fastqtl")
miscdir = file.path(datadir, "misc")
genodir = file.path(datadir, "geno")


## set some useful filenames
samps_file   = file.path(miscdir,"samps.txt")
deliv_file   = file.path(miscdir, "delivery_date.txt")
sex_age_file = file.path(miscdir, "covariates.sex_age")
peer_file    = file.path(rasqdir, "Sample_ALL.secondpass.gene.norm.tpm.10.txt")
pca_file     = file.path(genodir, "5M.imputed.dose.noambig.ldprune.hwe.extract.evec")
ctype_file   = file.path(miscdir, "all.secondpass.rsem.genes.noribo.LM22.cibersort.csv")
out_file     = file.path(fastdir, "Sample_ALL.fastqtl.covars.txt")


## get vector of samples to be analyzed
samps = readLines(samps_file)
samp_samps = paste("Sample", samps, sep="_")


## read in additional covariate informamtion
deliv_date   = read.delim(deliv_file, sep="\t", header=F, col.names=c("Sample", "date"), stringsAsFactors=F); row.names(deliv_date) = deliv_date$Sample
sex_age      = read.delim(sex_age_file, sep="", header=F, col.names=c("Samp","Sample","foo","sex","age"), row.names=2, stringsAsFactors=F)[samps,c("sex","age")]
peer_factors = read.delim(peer_file, sep="\t", header=F, row.names=1, col.names=paste("peer", 1:10, sep=""))[samps,]
pca_evec     = read.delim(pca_file, sep="", header=F, skip=1, row.names=1, stringsAsFactors=F, na.strings="???")
ctype_frac   = read.csv(ctype_file, header=T, row.names=1, stringsAsFactors=F)[samps,1:22]

ctype_frac.clean = ctype_frac[,apply(ctype_frac, 2, sum)>0]

## change row names because of odd formatting
row.names(pca_evec) = strsplit(row.names(pca_evec),":") %>% unlist %>% .[c(TRUE,FALSE)]
pca_evec = pca_evec[samps,]
colnames(pca_evec) = c(paste("pc", 1:10, sep=""), "foo")
         

## create covariate table
covars = t(cbind(sex_age, model.matrix( ~ date, deliv_date[samps,])[,-1], peer_factors[samps,], pca_evec[samps,1:3], ctype_frac.clean[samps,]))

## save covariates
write.table(covars, out_file, quote=F, sep="\t")
