#!user/bin/env Rscript

library(dplyr)
library(tidyr)



args = commandArgs(trailingOnly=TRUE)
pop  = args[1]

## set some useful directories
datadir = "/project/tishkofflab/rna_redo/afr_wbrna/data"
sampdir = file.path(datadir, "samples")
rasqdir = file.path(datadir, "rasqual")
fastdir = file.path(datadir, "fastqtl")
miscdir = file.path(datadir, "misc")
genodir = file.path(datadir, "geno")
gpopdir = file.path(genodir, pop)
fpopdir = file.path(fastdir, pop)
rpopdir = file.path(rasqdir, pop)


## set some useful filenames
samps_file   = paste(c(miscdir,"/pops/in5M/",pop,".txt"), collapse="")
deliv_file   = file.path(miscdir, "delivery_date.txt")
sex_age_file = file.path(miscdir, "covariates.sex_age")
peer_file    = file.path(rpopdir, paste(pop, "secondpass.gene.norm.tpm.10.txt", sep="."))
pca_file     = file.path(gpopdir, paste(pop, "5M.imputed.dose.noambig.ldprune.hwe.extract.evec", sep="."))
ctype_file   = file.path(miscdir, "all.secondpass.rsem.genes.noribo.LM22.cibersort.csv")
out_file     = file.path(fastdir, pop, paste(pop, "fastqtl.covars.txt", sep="."))


## get vector of samples to be analyzed
samps = readLines(samps_file)
samp_samps = paste("Sample", samps, sep="_")


## read in additional covariate informamtion
deliv_date   = read.delim(deliv_file, sep="\t", header=F, col.names=c("Sample", "date"), stringsAsFactors=F); row.names(deliv_date) = deliv_date$Sample
sex_age      = read.delim(sex_age_file, sep="", header=F, col.names=c("Samp","Sample","foo","sex","age"), row.names=2, stringsAsFactors=F)[samps,c("sex","age")]
peer_factors = read.delim(peer_file, sep="\t", header=F, row.names=1, col.names=c("sample",paste("peer", 1:10, sep="")))[samps,]
pca_evec     = read.delim(pca_file, sep="", header=F, skip=1, row.names=1, stringsAsFactors=F, na.strings="???")
ctype_frac   = read.csv(ctype_file, header=T, row.names=1, stringsAsFactors=F)[samps,1:22]

ctype_frac.clean = ctype_frac[,apply(ctype_frac, 2, sum)>0]

## change row names because of odd formatting
row.names(pca_evec) = strsplit(row.names(pca_evec),":") %>% unlist %>% .[c(TRUE,FALSE)]
pca_evec = pca_evec[samps,]
colnames(pca_evec) = c(paste("pc", 1:10, sep=""), "foo")
         

## create covariate table
covars = t(cbind(sex_age, model.matrix( ~ date, deliv_date[samps,])[,-1], peer_factors[samps,], pca_evec[samps,1:3], ctype_frac.clean[samps,]))
covars.df = cbind(data.frame(id=row.names(covars)), covars %>% as.data.frame)

## save covariates
write.table(covars.df, out_file, quote=F, row.names=F, sep="\t")
