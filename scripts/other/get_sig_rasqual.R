#!/usr/bin/env Rscript

library(data.table)
library(dplyr)
library(tidyr)



args = commandArgs(trailingOnly=TRUE)
rasqfile = args[1]
fdr_thresh = args[2]

## set some directories
datadir = "/project/tishkofflab/rna_redo/afr_wbrna/data"
rasqdir = file.path(datadir, "rasqual")
miscdir = file.path(datadir, "misc")

## useful data to have on hand
rasq_cols = readLines(file.path(miscdir, "rasqual_cols.txt"))


## read in rasqual data
rasq = na.omit(fread(rasqfile,
             col.names = rasq_cols,
             stringsAsFactors = F))
rasq[, p.value := pchisq(chisq, df=1, lower.tail=F),]
rasq[, bhp := p.adjust(p.value, method="BH"), by = gene_id]


## get FDR threshold: minimum BH-corrected p-value per gene. For ties, variants are randomly selected
rasq_bh_bh_thresh = rasq[, min(bhp), by = gene_id] %>%
    mutate(bh_bhp = p.adjust(V1, method="BH")) %>%
    filter(bh_bhp < 0.05) %>%
    select(V1) %>% max


## get eQTL with corrected p-value below the FDR threshold
rasq_thresh = rasq %>%
    filter(bhp <= rasq_bh_bh_thresh)


## write out results
write.table(rasq_thresh, "", row.names=F, sep="\t", quote=F)
