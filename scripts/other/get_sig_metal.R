#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)



args = commandArgs(trailingOnly=TRUE)
fdr_thresh = args[1]

## set some directories
datadir = "/project/tishkofflab/rna_redo/afr_wbrna/data"
miscdir = file.path(datadir, "misc")
metal_file = file.path(datadir, "metal/METAANALYSIS1.TBL")


## read in metalual data
metal = read.delim(file = metal_file,
                  stringsAsFactors = F) %>%
    filter(complete.cases(.)) %>%
    separate(MarkerName, into=c("gene_id", "chr", "pos"), sep=":", remove=F) %>%
    group_by(gene_id) %>%
    mutate(bhp = p.adjust(P.value, method="BH")) %>%
    as.data.frame


## get minimum BH-corrected p-value per gene. For ties, variants are randomly selected
metal_minbh = metal %>%
    group_by(gene_id) %>%
    summarise(min_bhp = min(bhp)) %>%
    group_by(gene_id) %>% 
    sample_n(1) %>%
    as.data.frame


## get the FDR threshold
metal_bh_bh_thresh = metal_minbh %>%
    mutate(bh_bhp = p.adjust(min_bhp, method="BH")) %>%
    filter(bh_bhp < fdr_thresh) %>%
    select(min_bhp) %>% max


## get eQTL with corrected p-value below the FDR threshold
metal_thresh = metal %>%
    filter(bhp <= metal_bh_bh_thresh)


## write out results
write.table(metal_thresh, "", row.names=F, sep="\t", quote=F)
