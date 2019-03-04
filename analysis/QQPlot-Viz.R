library(data.table)
library(tidyverse)
library(Homo.sapiens)
library(qqman)
library(qvalue)

setwd("/scratch/groups/rmccoy22/aseyedi2/QTLtoolsResults-WHLBLD")

dt <- fread("nominals.all.chunks.NE_only.txt.gz") %>%
  setnames(., c("intron_cluster", "chrom", "pheno_start", "pheno_end", 
                "strand", "total_cis", "distance", "variant_id", "variant_chrom", 
                "var_start", "var_end", "p", "beta", "is_top_variant", "drop")) %>%
  setorder(., p)

dt[, drop := NULL]

head(dt[is_top_variant == TRUE])

gene_list <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene) %>%
  as.data.table()
gene_list[, subjectHits := .I]
symbols <- org.Hs.egSYMBOL %>%
  as.data.table()
gene_list <- gene_list[symbols, on = "gene_id", nomatch = 0]
gene_list <- gene_list[, c("subjectHits", "symbol")]

coords_gr <- dt[, c("chrom", "pheno_start", "pheno_end")] %>%
  mutate(., chrom = paste0("chr", chrom)) %>%
  makeGRangesFromDataFrame

olaps <- findOverlaps(coords_gr, genes(TxDb.Hsapiens.UCSC.hg19.knownGene)) %>%
  as.data.table()

dt[, queryHits := .I]

dt <- dt[olaps, on = "queryHits", nomatch = 0]

dt <- dt[gene_list, on = "subjectHits", nomatch = 0]

setorder(dt, p)

qq(dt$p)

#########

gtp <- fread("permutations_full.txt.gz") %>%
  setnames(., c("intron_cluster", "chrom", "pheno_start", "pheno_end", 
                "strand", "total_cis", "distance", "variant_id", "variant_chrom", 
                "var_start", "var_end", "df", "dummy", "param_1", "param_2",
                "p", "beta", "emp_p", "adj_p")) %>%
  setorder(., adj_p)

gtp[, qval := qvalue(gtp$adj_p)$qvalue]

neand <- fread("/scratch/groups/rmccoy22/Ne_sQTL/sra/juncfiles/intronclustering/WHLBLD/tag_snps.neand.EUR.bed") %>%
  mutate(., var_id_1 = paste(V1, V3, V4, V5, "b37", sep = "_")) %>%
  mutate(., var_id_2 = paste(V1, V3, V5, V4, "b37", sep = "_")) %>%
  as.data.table()

neand_list <- c(neand$var_id_1, neand$var_id_2)

gtp[, is_neand := variant_id %in% neand_list]

gtp[qval < 0.1 & is_neand == TRUE]

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