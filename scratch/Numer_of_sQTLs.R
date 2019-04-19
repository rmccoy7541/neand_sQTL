library(data.table)
library(tidyverse)
library(qvalue)
library(pbapply)

setwd("~/Desktop/sQTL/")

count_sqtl <- function(tissue) {
  gtp <- fread(paste0(tissue, ".permutations_full.txt.gz")) %>%
    setnames(., c("intron_cluster", "chrom", "pheno_start", "pheno_end", 
                  "strand", "total_cis", "distance", "variant_id", "variant_chrom", 
                  "var_start", "var_end", "df", "dummy", "param_1", "param_2",
                  "p", "beta", "emp_p", "adj_p")) %>%
    setorder(., adj_p)
  
  gtp[, qval := qvalue(gtp$adj_p)$qvalues]
  
  neand <- fread("sprime_calls.txt.gz")[vindija_match == "match" | altai_match == "match"] %>%
    mutate(., var_id = paste(CHROM, POS, REF, ALT, "b37", sep = "_")) %>%
    as.data.table()
  
  gtp[, is_neand := variant_id %in% neand$var_id]
  
  return(data.table(TISSUE = tissue, n_sqtl = nrow(gtp[is_neand == TRUE & qval < 0.1])))
}

sqtl_counts_by_tissue <- do.call(rbind, lapply(c("BRNCHA", "BRNCTXB", "LIVER", "MSCLSK", "TESTIS", "WHLBLD"), function(x) count_sqtl(x)))
sqtl_counts_by_tissue[, n_samples := c(154, 118, 153, 491, 225, 369)]

ggplot(data = sqtl_counts_by_tissue, aes(x = n_samples, y = n_sqtl, label = TISSUE, color = TISSUE)) +
  theme_bw() +
  geom_point() +
  geom_text(hjust = 0, nudge_x = 5, nudge_y = 1) +
  xlim(0, 600) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  xlab("Number of Genotyped Samples") +
  ylab("Number of sQTL (FDR = 0.1)")