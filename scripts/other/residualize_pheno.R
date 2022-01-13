#!/usr/bin/Rscript

library(dplyr)

## read in args
args = commandArgs(trailingOnly=TRUE)

pheno_file = args[1]
covar_dir = args[2]
covar_pattern = args[3]
out_pref = args[4]

pheno = read.delim(pheno_file, stringsAsFactors=F)
samps = pheno %>% select(-c(X.Chr,start,end,ID)) %>% colnames
Y.mat = (pheno %>% select(samps) %>% t)

## list files
covar_files = list.files(path=covar_dir,
                         pattern=covar_pattern)

for (covar_file in covar_files){

    ## read in covariates for given number of peer factors
    covars = read.delim(file.path(covar_dir, covar_file), stringsAsFactors=F) %>% t

    ## get residuals by normal multiple regression
    residuals = Y.mat - covars %*% solve(t(covars) %*% covars) %*% t(covars) %*% Y.mat

    ## write results
    out.f = paste(out_pref, covar_file, sep="_")
    write.table(cbind(pheno %>% select(c(X.Chr,start,end,ID)), t(residuals)),
                file = out.f,
                sep = "\t",
                quote = F,
                row.names = F)
}
