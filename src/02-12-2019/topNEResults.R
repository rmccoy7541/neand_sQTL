library(data.table)
library(tidyverse)
library(qqman)

setwd("/scratch/groups/rmccoy22/Ne_sQTL/sra/juncfiles/intronclustering/WHLBLD")

dt <- do.call(rbind, lapply(1:100, function(x) tryCatch(fread(paste0("WHLBLD_nominals_chunk_", x, "_out.txt")), 
                                                        error = function(e) NULL))) %>%
  setnames(., c("intron_cluster", "chrom", "pheno_start", "pheno_end", 
                "strand", "total_cis", "distance", "variant_id", "variant_chrom", 
                "var_start", "var_end", "p", "beta", "is_top_variant", "drop")) %>%
  setorder(., p)

dt[, drop := NULL]

dt[is_top_variant == TRUE]

qq(dt$p)
