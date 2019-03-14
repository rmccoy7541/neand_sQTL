library(data.table)
library(tidyverse)
library(Homo.sapiens)
library(qqman)
library(qvalue)

install.packages("Homo.sapiens", "qqman", "qvalue")

setwd("/home/arta/Documents/GitHub/Ne-sQTL/data/QTLTools_Results/")


# installing/loading the latest installr package:
install.packages("installr"); library(installr) # install+load installr

updateR() # updating R.

gene_list <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene) %>%
  as.data.table()
gene_list[, subjectHits := .I]
symbols <- org.Hs.egSYMBOL %>%
  as.data.table()
gene_list <- gene_list[symbols, on = "gene_id", nomatch = 0]
gene_list <- gene_list[, c("subjectHits", "symbol")]

read_tissue_results <- function(tissue_id) {
  dt <- fread(paste0(tissue_id, ".permutations_full.txt.gz")) %>%
    setnames(., c("intron_cluster", "chrom", "pheno_start", "pheno_end", 
                  "strand", "total_cis", "distance", "variant_id", "variant_chrom", 
                  "var_start", "var_end", "df", "dummy", "param_1", "param_2",
                  "p", "beta", "emp_p", "adj_p")) %>%
    setorder(., adj_p)
  dt[, tissue := tissue_id]
  dt[, qval := qvalue(dt$adj_p)$qvalues]
  return(dt)
}

gtp <- do.call(rbind, lapply(c("BRNCHA", "TESTIS", "WHLBLD"), function(x) read_tissue_results(x)))

neand <- fread("tag_snps/tag_snps.neand.EUR.bed") %>%
  mutate(., var_id_1 = paste(V1, V3, V4, V5, "b37", sep = "_")) %>%
  mutate(., var_id_2 = paste(V1, V3, V5, V4, "b37", sep = "_")) %>%
  as.data.table()

neand_list <- c(neand$var_id_1, neand$var_id_2)

gtp[, is_neand := variant_id %in% neand_list]

nrow(gtp[qval < 0.1 & is_neand == TRUE])

coords_gr <- gtp[, c("chrom", "pheno_start", "pheno_end")] %>%
  mutate(., chrom = paste0("chr", chrom)) %>%
  makeGRangesFromDataFrame

olaps <- findOverlaps(coords_gr, genes(TxDb.Hsapiens.UCSC.hg19.knownGene)) %>%
  as.data.table()

gtp[, queryHits := .I]

gtp <- gtp[olaps, on = "queryHits", nomatch = 0]

gtp <- gtp[gene_list, on = "subjectHits", nomatch = 0]

setorder(gtp, adj_p)
gtp[qval < 0.1 & is_neand == TRUE]

setorder(gtp, symbol)
