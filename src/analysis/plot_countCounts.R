library(data.table)
library(tidyverse)
library(ggplot2)

setwd(snakemake@input[["countCountsDir"]])

file_list <- list.files(pattern = "_countCounts.txt", full.names = TRUE)
f_dowle2 = function(DT) {
  for (i in names(DT))
    DT[is.na(get(i)), (i):=0]
}
read_tissue_counts <- function(filename) {
  dt_tissue <- fread(filename, header = TRUE) 
  f_dowle2(dt_tissue)
  tissue_name <- gsub("^([^_]*_[^_]*)_.*$", "\\1", filename)
  tissue_name <- gsub("/Users/rajivmccoy/Downloads/", "", tissue_name)
  dt_tissue[, tissue := tissue_name]
  return(dt_tissue[!duplicated(transcript_id)])
}
dt <- do.call(rbind, lapply(file_list, function(x) read_tissue_counts(x)))

png(filename = "results/SeparateTissues.png")
ggplot(dt, aes(x = `HH>0`, y = `HN>0` + `NN>0`)) +
  geom_bin2d() + 
  scale_fill_gradient(name = "count", trans = "log10") +
  facet_wrap(~ tissue)
def.off()

png(filename = "results/AllTissuesNLIso.png")
ggplot(dt, aes(x = `HH>0`, y = `HN>0` + `NN>0`)) +
  geom_bin2d(bins = 100) + 
  scale_fill_gradient(name = "count", trans = "log10")
def.off()