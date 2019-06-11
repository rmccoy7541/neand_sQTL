library(data.table)
library(tidyverse)
library(qvalue)
library(pbapply)
library(Homo.sapiens)


args = commandArgs(trailingOnly=TRUE)
# 1 is WD, 2 is sprime file


# Working directory as commandline arg, make sure all of the catt'd permutation files are in your dir
setwd(args[1])

count_sqtl <- function(tissue, summarize = FALSE) {
  gtp <- fread(paste0(tissue, "_permutations.txt")) %>%
    setnames(., c("intron_cluster", "chrom", "pheno_start", "pheno_end", 
                  "strand", "total_cis", "distance", "variant_id", "variant_chrom", 
                  "var_start", "var_end", "df", "dummy", "param_1", "param_2",
                  "p", "beta", "emp_p", "adj_p")) %>%
    setorder(., adj_p)
  
  gtp[, qval := qvalue(gtp$adj_p)$qvalues]
  

  neand <- fread(args[2])[vindija_match == "match" | altai_match == "match"] %>%
    mutate(., var_id = paste(CHROM, POS, REF, ALT, "b37", sep = "_")) %>%
    as.data.table()
  
  gtp[, logP := -log10(adj_p)]
  setorder(gtp, logP)
  gtp[, expectedP := rev(-log10(ppoints(n = length(gtp$adj_p))))]
  
  gtp[, is_neand := variant_id %in% neand$var_id]
  
  if (summarize == TRUE) {
    return(data.table(TISSUE = tissue, n_sqtl = nrow(gtp[is_neand == TRUE & qval < 0.1]), n_total = nrow(gtp[qval < 0.1])))
  } else {
    return(data.table(gtp, TISSUE_ID = tissue))
  }
  
}

sqtl_counts_by_tissue <- do.call(rbind, lapply(c("ADPSBQ","ADPVSC","ADRNLG","ARTACRN","ARTAORT","ARTTBL","BREAST","BRNACC","BRNAMY","BRNCDT","BRNCHA","BRNCHB","BRNCTXA","BRNCTXB","BRNHPP","BRNHPT","BRNNCC","BRNPTM","BRNSNG","BRNSPC","CLNSGM","CLNTRN","ESPGEJ","ESPMCS","ESPMSL","FIBRLBLS","HRTAA","HRTLV","LCL","LIVER","LUNG","MSCLSK","NERVET","OVARY","PNCREAS","PRSTTE","PTTARY","SKINNS","SKINS","SLVRYG","SNTTRM","SPLEEN","STMACH","TESTIS","THYROID","UTERUS","VAGINA","WHLBLD"), function(x) count_sqtl(x, summarize = TRUE)))
sqtl_counts_by_tissue[, n_samples := c(385,313,175,152,267,388,251,109,88,144,154,125,136,118,111,108,130,111,80,83,203,246,213,358,335,300,264,272,117,153,383,491,361,122,220,132,157,335,414,85,122,146,237,225,399,101,106,369)]

ggplot(data = sqtl_counts_by_tissue, aes(x = n_samples, y = n_sqtl, label = TISSUE, color = TISSUE)) +
  theme_bw() +
  geom_point() +
  ggrepel::geom_text_repel() + # hjust = 0, nudge_x = 4, nudge_y = 2
  xlim(0, 600) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  xlab("Number of Genotyped Samples") +
  ylab("Number of sQTL (FDR = 0.1)")
  
tissue_gtp <- do.call(rbind, lapply(c("ADPSBQ","ADPVSC","ADRNLG","ARTACRN","ARTAORT","ARTTBL","BREAST","BRNACC","BRNAMY","BRNCDT","BRNCHA","BRNCHB","BRNCTXA","BRNCTXB","BRNHPP","BRNHPT","BRNNCC","BRNPTM","BRNSNG","BRNSPC","CLNSGM","CLNTRN","ESPGEJ","ESPMCS","ESPMSL","FIBRLBLS","HRTAA","HRTLV","LCL","LIVER","LUNG","MSCLSK","NERVET","OVARY","PNCREAS","PRSTTE","PTTARY","SKINNS","SKINS","SLVRYG","SNTTRM","SPLEEN","STMACH","TESTIS","THYROID","UTERUS","VAGINA","WHLBLD"), function(x) count_sqtl(x)))

ggplot(data = tissue_gtp[is_neand == TRUE], aes(x = expectedP, y = logP, color = TISSUE_ID)) +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.position = "none") +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) +
  xlab(expression(Expected -log[10](italic("p")))) +
  ylab(expression(Observed -log[10](italic("p")))) +
  facet_wrap(~ TISSUE_ID)

gene_list <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene) %>%
  as.data.table()
gene_list[, subjectHits := .I]
symbols <- org.Hs.egSYMBOL %>%
  as.data.table()
gene_list <- gene_list[symbols, on = "gene_id", nomatch = 0]
gene_list <- gene_list[, c("subjectHits", "symbol")]

coords_gr <- tissue_gtp[, c("chrom", "pheno_start", "pheno_end")] %>%
  mutate(., chrom = paste0("chr", chrom)) %>%
  makeGRangesFromDataFrame

olaps <- findOverlaps(coords_gr, genes(TxDb.Hsapiens.UCSC.hg19.knownGene)) %>%
  as.data.table()

tissue_gtp[, queryHits := .I]

tissue_gtp_annotated <- tissue_gtp[olaps, on = "queryHits"]

tissue_gtp_annotated <- tissue_gtp_annotated[gene_list, on = "subjectHits"] %>%
  setorder(., adj_p)

dplyr::select(head(tissue_gtp_annotated[qval < 0.1 & is_neand == TRUE], 100), intron_cluster, variant_id, TISSUE_ID, adj_p, qval, symbol)