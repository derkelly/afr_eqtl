#!/usr/bin/Rscript

args = commandArgs(trailingOnly=TRUE)
input = args[1]

input = file('stdin', 'r')

values = read.delim(input, sep="\t", header=F. stringsAsFactors=F)

print(apply(values,1,sum))
