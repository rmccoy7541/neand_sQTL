library(data.table)
library(tidyverse)
library(qqman)
library(qvalue)
library(org.Hs.eg.db)
library(annotate)
library(rtracklayer)
# commVars[1] is GTEx perm pass result, "gtf" is gencode.v31.annotation.gtf, "sprime" is sprime_calls.txt

#########
gene_list <- rtracklayer::import(snakemake@input["gtf"]) %>%
  makeGRangesFromDataFrame(., keep.extra.columns = T) %>%
  as.data.table()

gene_list <- gene_list[type == "gene" & gene_type == "protein_coding"]

gene_list[, subjectHits := .I]

# Perm Pass QTLTOOLS
gtp <- fread(snakemake@input["perm"]) %>%
  setorder(., pval_nominal)

# Not TAGSNPS but SPRIME
neand <- fread(snakemake@input["sprime"])[vindija_match == "match" | altai_match == "match"] %>%
  mutate(., var_id_1 = paste(CHROM, POS, REF, ALT, "b38", sep = "_")) %>%
  as.data.table()

neand_list <- c(neand$var_id_1, neand$var_id_2)

gtp[, is_neand := variant_id %in% neand_list]

gtp <- gtp[is_neand == TRUE]

coords_tab <- strsplit(gtp$phenotype_id, split = ":")

coords_gr <- as.data.frame(do.call(rbind, coords_tab))

names(coords_gr) <- c("chr", "start", "end")

coords_gr <- makeGRangesFromDataFrame(coords_gr, keep.extra.columns = F)
gene_list_gr <- makeGRangesFromDataFrame(gene_list, keep.extra.columns = T)

olaps <- findOverlaps(coords_gr, gene_list_gr) %>%
  as.data.table()

gtp[, queryHits := .I]

gtp <- gtp[olaps, on = "queryHits", nomatch = 0]

gtp <- gtp[gene_list, on = "subjectHits", nomatch = 0]

setorder(gtp, pval_nominal)
table <- gtp[is_neand == TRUE]

tis_name <- strsplit("Whole_Blood.v8.sqtl_signifpairs.txt", split = "[.]")[[1]][1]

write.table(table, paste0(tis_name, "_permutation_table_NE.txt"), row.names=F, quote=F, sep="\t")

