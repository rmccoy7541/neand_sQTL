library(data.table)
library(GenomicRanges)
library(stringr)
library(biomaRt)
library(annotate)
library(rtracklayer)
library(dplyr)

# use this file which contains information about all of the sequences in hg38
gene_list <- rtracklayer::import("C:/Users/artas/Downloads/gencode.v26.GRCh38.genes.gtf") %>%
  makeGRangesFromDataFrame(., keep.extra.columns = T) %>%
  as.data.table()

# subset to just protein coding genes
gene_list <- gene_list[type == "gene" & gene_type == "protein_coding"]

gene_list[, subjectHits := .I]

# read coords for NL-specificish introns
dt <- fread("results/loosenedRestrictions.txt", header = T)

coords <- as.data.table(str_split_fixed(dt$transcript_id, "_", 3))

names(coords) <- c("seqnames", "start", "stop")
  
coords_gr <- makeGRangesFromDataFrame(coords, keep.extra.columns = F)
gene_list_gr <- makeGRangesFromDataFrame(gene_list, keep.extra.columns = T)

# find overlaps like in the top genes script
olaps <- findOverlaps(coords_gr, gene_list_gr) %>%
  as.data.table()

dt[, queryHits := .I]

dt <- dt[olaps, on = "queryHits", nomatch = 0]

dt <- dt[gene_list, on = "subjectHits", nomatch = 0]

dt <- select(dt, c(transcript_id, seqnames, start, end, gene_name, n))

write.table(dt,
            file = paste0("loosenedRestrictionsGenes.txt"),
            sep = "\t",
            row.names = F,
            quote = FALSE)