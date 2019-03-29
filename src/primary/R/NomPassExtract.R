library(data.table)
library(tidyverse)
commvars = commandArgs(trailingOnly=TRUE)
# basedir <- "/scratch/groups/rmccoy22/Ne_sQTL/sra/juncfiles/intronclustering/WHLBLD/"

# Whole blood chunk 1-100 is arg 1
sqtl <- fread(commvars[1]) %>% 
  setnames(., c("intron_cluster", "chrom", "pheno_start", "pheno_end", "strand", "total_cis", "distance", "variant_id", "variant_chrom", "var_start", "var_end", "p", "beta", "is_top_variant"))

# Neanderthal bed file is arg 2
neand <- fread(commvars[2]) %>%
  mutate(., var_id_1 = paste(V1, V3, V4, V5, "b37", sep = "_")) %>%
  mutate(., var_id_2 = paste(V1, V3, V5, V4, "b37", sep = "_")) %>%
  as.data.table()
  
neand_list <- c(neand$var_id_1, neand$var_id_2)

sqtl[, is_neand := variant_id %in% neand_list]

neand_sqtl <- sqtl[is_neand == TRUE] %>%
  setorder(., p)

name <- tools::file_path_sans_ext(commvars[1])

write.table(neand_sqtl, file = paste0(name, "_out.txt"), quote = F, col.names = F, row.names = F)
