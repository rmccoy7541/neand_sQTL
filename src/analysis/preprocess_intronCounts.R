library(tidyr)
library(dplyr)
library(data.table)

introns <-
  fread(
    snakemake@input["introns"],
    stringsAsFactors = FALSE,
    header = TRUE
  )

tistab <- fread(snakemake@input["tistab"])

sites <- with(tistab, split(SAMPID, TISSUE))

keep <- c('Name', 'Description')

tissues <- lapply(sites, function(x)
  introns[, .SD, .SDcols = c(keep, intersect(names(introns), x))])

sapply(names(tissues_kept), function (x) write.table(tissues_kept[[x]], file=paste0("results/IC/",x,"_intronCounts.txt"), row.names=F, quote=FALSE, sep="\t"))
