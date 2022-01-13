#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)


args = commandArgs(trailingOnly=TRUE)
sqtlfile = args[1]
fdr_thresh = args[2]


## set some directories
datadir = "/project/tishkofflab/rna_redo/afr_wbrna/data"
leafdir = file.path(datadir, "leafcutter")
resultdir = file.path(datadir, "leafcutter")
miscdir = file.path(datadir, "misc")

## useful data to have on hand
sqtl_cols = readLines(file.path(miscdir, "gemma_cols.lmm4.txt"))


## read in sqtl data
sqtl = read.delim(file = sqtlfile,
                  col.names = sqtl_cols,
                  stringsAsFactors = F,
                  na.strings = "?") %>%
    filter(complete.cases(.)) %>%
    group_by(gene_id) %>%
    mutate(bhp = p.adjust(p.value, method="BH")) %>%
    as.data.frame


## get minimum BH-corrected p-value per gene. For ties, variants are randomly selected
sqtl_minbh = sqtl %>%
    group_by(gene_id) %>%
    summarise(min_bhp = min(bhp)) %>%
    group_by(gene_id) %>% 
    sample_n(1) %>%
    as.data.frame


## get the FDR threshold
sqtl_bh_bh_thresh = sqtl_minbh %>%
    mutate(bh_bhp = p.adjust(min_bhp, method="BH")) %>%
    filter(bh_bhp < fdr_thresh) %>%
    select(min_bhp) %>% max


## get eQTL with corrected p-value below the FDR threshold
sqtl_thresh = sqtl %>%
    filter(bhp <= sqtl_bh_bh_thresh)


## write out results
write.table(sqtl_thresh, "", row.names=F, sep="\t", quote=F)
