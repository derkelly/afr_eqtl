#!/usr/bin/Rscript

source("general_fncs.R")

library(dplyr)
library(tidyr)
library(peer)

leafdir = "/project/tishkofflab/rna_redo/afr_wbrna/data/leafcutter"
fastdir = "/project/tishkofflab/rna_redo/afr_wbrna/data/fastqtl"
miscdir = "/project/tishkofflab/rna_redo/afr_wbrna/data/misc"


## read in normalized data
qqnorm = read.delim(file.path(leafdir, "pheno", "leafcutter_all_perind.chr.counts.qqnorm_all.geno.gz"),
stringsAsFactors=F)


## read in other covariates, previously used by fastQTL
covars = read.delim(file.path(fastdir, "Sample_ALL.fastqtl.covars.txt"), stringsAsFactors=F, row.names=1)


## get peer factors, 5, 10, 15, 20, 25, and 30.
set.seed(4690)

n.factors.vec = c(1:20)

for (n.factors in n.factors.vec){

qqnorm.peer = run_peer(t(qqnorm %>% select(-c("X.Chr","start","end","ID"))), num_factors=n.factors)
colnames(qqnorm.peer) = paste("peer", 1:n.factors, sep="")
row.names(qqnorm.peer) = qqnorm %>% select(-c("X.Chr","start","end","ID")) %>% colnames


## remove previous peer factors, add new peer factors, and write to covariate file.
covars.t = covars %>%
   t %>%
   as.data.frame %>%
   select(-starts_with("peer"))


# ## format leafcutter fractions as a bed file
# leaf.frac.clean.bed = leaf.frac.clean %>% 
# as.data.frame %>% 
# rownames_to_column(., var="ID") %>% 
# separate(ID, into=c("#Chr","start","end","cluster"), sep=":", remove=F) %>%
# select(c("#Chr","start","end","ID",row.names(covars.t)))


## add new peer factors, making sure that the columns match
covars.t[,paste("peer", 1:n.factors, sep="")] = qqnorm.peer[row.names(covars.t),]
write.table(t(covars.t), file=file.path(leafdir, paste(c("Sample_ALL.leafcutter.covars.peer",n.factors,".txt"), collapse="")), quote=F, row.names=F, sep="\t")
}
