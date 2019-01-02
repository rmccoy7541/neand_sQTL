#!/usr/bin/env Rscript
require("data.table")
# arg 1 is leafcutter PCs, arg 2 is GTEx covariate file for tissue
args = commandArgs(trailingOnly=TRUE)
library("data.table")
leafcut <- fread(args[1])
#leafcut <- fread("testNE_sQTL_perind.counts.gz.PCs")
gtexPC <- fread(args[2])
#gtexPC <- fread("GTEx_Analysis_v7_eQTL_covariates/Whole_Blood.v7.covariates.txt")
# make the use of "id" header standard
setnames(gtexPC, "ID", "id")
# concatenate; fill missing values with na, but then remove 
out = rbind(leafcut,gtexPC, use.names=T, fill=T)
out <- t(na.omit(t(out)))
tissue = tools::file_path_sans_ext(basename(args[2]))
outfile <- gsub(" ", "", paste(tissue, "_output.txt"), fixed = TRUE)
write.table(out, file = outfile, quote = F, sep= "\t", eol = "\r\n", row.names = F)
