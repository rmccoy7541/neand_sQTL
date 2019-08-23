library(data.table)
library(tidyverse)
# commvars = commandArgs(trailingOnly=TRUE)

# Whole blood chunk 1-100 is arg 1
sqtl <- fread(snakemake@input[[1]]) %>% 
  setnames(., c("intron_cluster", "chrom", "pheno_start", "pheno_end", "strand", "total_cis", "distance", "variant_id", "variant_chrom", "var_start", "var_end", "p", "beta", "is_top_variant"))

# Neanderthal bed file is arg 2
neand <- fread(snakemake@input[[2]]) %>%
  mutate(., var_id_1 = paste(CHROM, POS, REF, ALT, "b37", sep = "_")) %>%
  as.data.table()
  
neand_list <- c(neand$var_id_1, neand$var_id_2)

sqtl[, is_neand := variant_id %in% neand_list]

neand_sqtl <- sqtl[is_neand == TRUE] %>%
  setorder(., p)

name <- tools::file_path_sans_ext(snakemake@input[[1]])

write.table(neand_sqtl, file = paste0(name, "_out.txt"), quote = F, col.names = F, row.names = F)
