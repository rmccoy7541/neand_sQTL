#!/usr/bin/env Rscript
require("data.table")
# arg 1 is leafcutter PCs, arg 2 is GTEx covariate file for tissue, 3 is tissue_table.txt, 4 is wd
args = commandArgs(trailingOnly=TRUE)
library("data.table")
leafcut <- fread(args[1], colClasses = "character")
#leafcut <- fread(as.character("Ne-sQTL_perind.counts.gz.PCs"))
gtexPC <- fread(args[2], colClasses = "character")
#gtexPC <- fread("GTEx_Analysis_v7_eQTL_covariates/Whole_Blood.v7.covariates.txt")
lookup <- fread(args[3], colClasses = "character")
# make the use of "id" header standard
setnames(gtexPC, "ID", "id")

# making sure our leafcutter friends have the right header names
ind <- match(names(leafcut), lookup$Run)
names(leafcut) <- lookup$submitted_subject_id[ind]

# concatenate; fill missing values with na, but then remove
out = rbind(leafcut,gtexPC, use.names=T, fill=T)
# omit NA's with a function that only works on rows
out <- t(na.omit(t(out)))
tissue = tools::file_path_sans_ext(basename(args[2]))
outfile <- gsub(" ", "", paste(tissue, "_output.txt"), fixed = TRUE)
write.table(out, file = outfile, quote = F, sep= "\t", eol = "\r\n", row.names = F)