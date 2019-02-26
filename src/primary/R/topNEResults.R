library(data.table)
library(tidyverse)
library(Homo.sapiens)
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
