#!user/bin/env Rscript

library(dplyr)
library(tidyr)



args = commandArgs(trailingOnly=TRUE)
pop = args[1]
outdir = args[2]


## set some useful directories
datadir = "/project/tishkofflab/rna_redo/afr_wbrna/data"
sampdir = paste(datadir, "samples", sep="/")
rasqdir = paste(datadir, "rasqual", sep="/")
miscdir = paste(datadir, "misc", sep="/")
genodir = paste(datadir, "geno", sep="/")
gpopdir = paste(genodir, pop, sep="/")
rpopdir = paste(rasqdir, pop, sep="/")


## set some useful filenames
samps_file   = paste(c(miscdir,"/pops/in5M/",pop,".txt"), collapse="")
deliv_file   = paste(miscdir, "delivery_date.txt", sep="/")
sex_age_file = paste(miscdir, "covariates.sex_age", sep="/")
peer_file    = paste(rpopdir, paste(pop, "secondpass.gene.norm.tpm.10.txt", sep="."), sep="/")
pca_file     = paste(gpopdir, paste(pop, "5M.imputed.dose.noambig.ldprune.hwe.extract.evec", sep="."), sep="/")
ctype_file   = paste(miscdir, "all.secondpass.rsem.genes.noribo.LM22.cibersort.csv", sep="/")


## get vector of samples to be analyzed
samps_all = readLines(paste(miscdir, "samples.txt", sep="/"))
samps = readLines(samps_file)
samp_samps = paste("Sample", samps, sep="_")


## read in additional covariate informamtion
deliv_date   = read.delim(deliv_file, sep="\t", header=F, col.names=c("Sample", "date"), stringsAsFactors=F); row.names(deliv_date) = deliv_date$Sample
sex_age      = read.delim(sex_age_file, sep="", header=F, col.names=c("Samp","Sample","foo","sex","age"), row.names=2, stringsAsFactors=F)[samps,c("sex","age")]
peer_factors = read.delim(peer_file, sep="\t", header=F, row.names=1)[samps,]
pca_evec     = read.delim(pca_file, sep="", header=F, skip=1, row.names=1, stringsAsFactors=F, na.strings="???")
ctype_frac   = read.csv(ctype_file, header=T, row.names=1, stringsAsFactors=F)[samps,1:22]

ctype_frac.clean = ctype_frac[,apply(ctype_frac, 2, sum)>0]

## change row names because of odd formatting
row.names(pca_evec) = strsplit(row.names(pca_evec),":") %>% unlist %>% .[c(TRUE,FALSE)]
pca_evec = pca_evec[samps,]

## create covariate table
covars = cbind(sex_age, model.matrix( ~ date, deliv_date[samps,])[,-1], peer_factors[samps,], pca_evec[samps,1:3], ctype_frac.clean[samps,])

## save covariates
out_file = paste(outdir, paste(pop, "X.txt", sep="."), sep="/")
write.table(covars, out_file, quote=F, sep="\t", row.names=F, col.names=F)
