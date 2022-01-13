#!/usr/bin/Rscript


args = commandArgs(trailingOnly=TRUE)
input = args[1]


fst = read.delim(input, sep="\t", header=T, na.strings="-nan")
fst[is.na(fst[,3]),3] = 0


print(cbind(fst[,c(1:2)],(fst[,3]-mean(fst[,3]))/(sd(fst[,3]))), row.names=F, col.names=F)
