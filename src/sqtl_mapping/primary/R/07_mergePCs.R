#!/usr/bin/env Rscript
require("data.table")
# arg 1 is leafcutter PCs, arg 2 is GTEx covariate file for tissue, 3 is tissue_table.txt
# args = commandArgs(trailingOnly=TRUE)
library("data.table")
leafcut <- fread(snakemake@input[[1]], colClasses = "character")
#leafcut <- fread(as.character("Ne-sQTL_perind.counts.gz.PCs"))
gtexPC <- fread(snakemake@input[[2]], colClasses = "character")
#gtexPC <- fread("Adipose_Subcutaneous.v7.covariates.txt")
lookup <- fread(snakemake@input[[3]], colClasses = "character")
#lookup <- fread("tissue_table.txt")
# make the use of "id" header standard
setnames(gtexPC, "ID", "id")

names(leafcut) <- sapply(strsplit(names(leafcut), "\\."), '[', 1)

# making sure our leafcutter friends have the right header names
ind <- match(names(leafcut), lookup$Sample_Name)
names(leafcut) <- lookup$submitted_subject_id[ind]

# concatenate; fill missing values with na, but then remove
out <- rbind(leafcut,gtexPC, use.names=T, fill=T)
# omit NA's with a function that only works on rows
out <- t(na.omit(t(out)))
tissue = tools::file_path_sans_ext(basename(args[2]))
#tissue <- tools::file_path_sans_ext(basename("Whole_Blood.v7.covariates.txt"))
outfile <- gsub(" ", "", paste(tissue, "_output.txt"), fixed = TRUE)
write.table(out, file = outfile, quote = F, sep= "\t", eol = "\r\n", row.names = F)
