#!/usr/bin/env Rscript

library(dplyr)
library(tidyr)



args = commandArgs(trailingOnly=TRUE)
pop = args[1]
outdir = args[2]


## set some useful directories
datadir = "/project/tishkofflab/rna_redo/afr_wbrna/data"
sampdir = paste(datadir,"samples",sep="/")
miscdir = paste(datadir,"misc",sep="/")

samps_file    = file.path(miscdir,"/pops/in5M/",paste(pop,"txt",sep="."))
gc_file       = file.path(miscdir, "secondpass.featureCounts.gc_content.txt")
feat_gex_file = file.path(sampdir,"Sample_ALL","Sample_ALL.secondpass.featurecounts.merge.txt")
tpm_gex_file  = file.path(sampdir,"Sample_ALL","Sample_ALL.secondpass.rsem.genes.results.tpm")

# get vector of samples to be analyzed
samps_all = readLines(file.path(miscdir, "samples.txt"))
samps = readLines(samps_file)
samp_samps = paste("Sample", samps, sep="_")

# read in the gc content information
gc_correct = read.delim(gc_file,sep="\t",header=T,stringsAsFactors=F)


# read in gene expression data and add gc content information
tpm      = read.delim(tpm_gex_file, sep="\t", header=T, stringsAsFactors=F)
features = read.delim(feat_gex_file, sep="\t", header=T, stringsAsFactors=F)
tpm.above_thresh = data.frame(gene_id=tpm$gene_id, above_thresh=apply(tpm[,samps], 1, function(x) sum(x > 0.1)))

features.clean = merge(features, tpm.above_thresh, by="gene_id")
features.clean = features.clean %>% filter(above_thresh > 4)
features.clean = merge(features.clean,gc_correct,by="gene_id")
features.clean$gc_frac = features.clean$gc_bp/features.clean$total_bp


# get bin bounds for gc correction
gc_bins = quantile(features.clean$gc_frac,seq(0,1,0.005))

# assign bins for gc values
features.clean$gc_bin = .bincode(features.clean$gc_frac,gc_bins)

# make long format for processing downstream
features.long = features.clean %>%
    gather(key=Sample, value=count, samps)

## add population column
#features.long = merge(features.long,samp_pop,by="Sample")


###############################
##### Relative Enrichment #####
###############################
#
# This is the calculation of the offset terms, K, for each individual
# in each population. We also need to calculate the GC correction term


## function to fit a spline to the data
fit_gc_spline = function(df){

    ## calculate relevant sums for GC correction for RASQUAL
    s_il = df %>%
        group_by(Sample,gc_bin) %>%
        summarise(s_il = sum(count))
    
    s_dotl = df %>%
        group_by(gc_bin) %>%
        summarise(s_dotl = sum(count))
    
    s_idot = df %>%
        group_by(Sample) %>%
        summarise(s_idot = sum(count))
    
    s_dotdot = sum(as.numeric(df$count))
    
    
    s_data = merge(s_il, s_dotl, by="gc_bin")
    s_data = merge(s_data, s_idot, by="Sample")
    
    
    s_data$f_il = log2( (s_data$s_il/s_data$s_dotl) / (s_data$s_idot/s_dotdot) )
    
    
    ## get mean gc_frac value for each bin
    mean_gc = df %>%
        group_by(gc_bin) %>%
        summarise(gc_frac_mean = mean(gc_frac))
    
    
    s_data = merge(s_data, mean_gc, by="gc_bin")
    s_data = s_data[complete.cases(s_data),]
    
    
    ## fit a smoothing spline
    fit.spline = smooth.spline(s_data$gc_frac_mean, s_data[complete.cases(s_data),]$f_il, spar=1)

    return(fit.spline)
}


spline_fits = fit_gc_spline(features.long)


## calculate normalized K and normalized feature counts
offsets = list()
norm_features = list()
norm_features.wide = list()


# get the number of features
N = length(unique(features.long$Sample))
J = length(unique(features.long$gene_id))

Y_dotdot = sum(as.numeric(features.long$count))/N/J

# calculate the offset term for each individual
Y_idot = features.long %>%
    group_by(Sample) %>%
    summarise(y_idot = sum(as.numeric(count))/J)

Y_idot$K = Y_idot$y_idot / Y_dotdot

features.long = merge(features.long, Y_idot[,c("Sample","K")], by="Sample")

features.long$K_gc = features.long$K * exp( spline_fits$y[features.long$gc_bin] )


# calculate offsets
offsets = spread(features.long[,c("gene_id","Sample","K_gc")], key=Sample, value=K_gc)



## write out data I'll need
K.fname = paste(outdir, paste(pop,"K.gc_correct.txt",sep="."), sep="/")
Y.fname = paste(outdir, paste(pop,"Y.txt",sep="."), sep="/")
    
write.table(offsets[,c("gene_id",samps)], K.fname, sep="\t", quote=F, row.names=F, col.names=F)
write.table(features.clean[,c("gene_id",samps)], Y.fname, sep="\t", quote=F, row.names=F, col.names=F)
