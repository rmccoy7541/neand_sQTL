library(data.table)
library(tidyverse)
library(qvalue)
library(pbapply)
library(Homo.sapiens)
library(pals)
library(qqman)
library(cowplot)
library(ggplotify)
library(ggrepel)
library(RColorBrewer)

count_sqtl <- function(tissue, summarize = FALSE, select_neand = TRUE) {
  gtp <- fread(paste0(snakemake@input[[1]]))
  
  neand <- fread(snakemake@input[[2]])[vindija_match == "match" | altai_match == "match"] %>%
    mutate(., var_id_1 = paste(CHROM, POS, REF, ALT, "b38", sep = "_")) %>%
    as.data.table()
  
  gtp[, is_neand := variant_id %in% neand$var_id]
  
  gtp_neand_t <- gtp[is_neand == TRUE]
  gtp_neand_f <- gtp[is_neand == FALSE]
  
  gtp_neand_t[, logP := -log10(adj_p)]
  setorder(gtp_neand_t, logP)
  gtp_neand_t[, expectedP := rev(-log10(ppoints(n = length(gtp_neand_t$adj_p))))]
  
  gtp_neand_f[, logP := -log10(adj_p)]
  setorder(gtp_neand_f, logP)
  gtp_neand_f[, expectedP := rev(-log10(ppoints(n = length(gtp_neand_f$adj_p))))]
  
  gtp <- rbind(gtp_neand_t, gtp_neand_f)
  
  if (summarize == TRUE) {
    return(data.table(TISSUE = tissue, n_sqtl = nrow(gtp[is_neand == select_neand & qval < 0.1]), n_total = nrow(gtp[qval < 0.1])))
  } else {
    return(data.table(gtp, TISSUE_ID = tissue))
  }
  
}

tissue_vector <- c("ADPSBQ","ADPVSC","ADRNLG","ARTACRN","ARTAORT","ARTTBL","BREAST","BRNACC","BRNAMY","BRNCDT","BRNCHA","BRNCHB","BRNCTXA","BRNCTXB","BRNHPP","BRNHPT","BRNNCC","BRNPTM","BRNSNG","BRNSPC","CLNSGM","CLNTRN","ESPGEJ","ESPMCS","ESPMSL","FIBRLBLS","HRTAA","HRTLV","LCL","LIVER","LUNG","MSCLSK","NERVET","OVARY","PNCREAS","PRSTTE","PTTARY","SKINNS","SKINS","SLVRYG","SNTTRM","SPLEEN","STMACH","TESTIS","THYROID","UTERUS","VAGINA","WHLBLD")

sqtl_counts_by_tissue <- do.call(rbind, lapply(tissue_names, function(x) count_sqtl(x, summarize=TRUE)))
sqtl_counts_by_tissue[, n_samples := c(763,564,275,450,253,770,177,213,291,263,298,325,425,243,236,277,232,182,164,480,527,192,389,432,401,622,559,452,689,100,251,867,181,1132,722,195,360,301,262,638,849,193,260,381,406,812,166,173,3288)]

tissue_table <- data.table(TISSUE_ID = tissue_vector, SAMPLE_SIZE = c(763,564,275,450,253,770,177,213,291,263,298,325,425,243,236,277,232,182,164,480,527,192,389,432,401,622,559,452,689,100,251,867,181,1132,722,195,360,301,262,638,849,193,260,381,406,812,166,173,3288))

sqtl_count_plot <- ggplot(data = sqtl_counts_by_tissue, aes(x = n_samples, y = n_sqtl, label = TISSUE, color = TISSUE)) +
  theme_bw() +
  geom_smooth(method = "lm", color = "darkgray", se = FALSE, lty = "dotted", fullrange = TRUE) +
  geom_point() +
  ggrepel::geom_text_repel() + # hjust = 0, nudge_x = 4, nudge_y = 2
  xlim(0, 600) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  xlab("Number of genotyped samples") +
  ylab("Number of sQTL (FDR = 0.1)") +
  scale_color_manual(values = colorRampPalette(brewer.pal(8, "Dark2"))(48)) +
  NULL

tissue_gtp <- do.call(rbind, lapply(tissue_vector, function(x) count_sqtl(x)))
tissue_gtp <- merge(tissue_gtp, tissue_table, "TISSUE_ID")

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

tissue_gtp_annotated <- merge(tissue_gtp, olaps, "queryHits", all.x = TRUE)

tissue_gtp_annotated <- merge(tissue_gtp_annotated, gene_list, "subjectHits", all.x = TRUE) %>%
  setorder(., adj_p)

to_file <- tissue_gtp_annotated[is_neand == TRUE][, -c("subjectHits", "queryHits")] %>%
  group_by(., tmp = paste(intron_cluster, chrom, pheno_start, pheno_end, strand, total_cis, distance, variant_id, TISSUE_ID)) %>%
  summarise(., phenotype_id = unique(intron_cluster), chrom = unique(chrom), pheno_start = unique(pheno_start), pheno_end = unique(pheno_end), 
            strand = unique(strand), total_cis = unique(total_cis), distance = unique(distance),
            variant_id = unique(variant_id), df = unique(df), dummy = unique(dummy), param_1 = unique(param_1), param_2 = unique(param_2), 
            p = unique(p), beta = unique(beta), emp_p = unique(emp_p), adj_p = unique(adj_p), qval = unique(qval), symbol = list(symbol)) %>%
  as.data.table() %>%
  setorder(., adj_p)

dplyr::select(head(tissue_gtp_annotated[qval < 0.1 & is_neand == TRUE], 100), intron_cluster, variant_id, TISSUE_ID, adj_p, qval, symbol)

tissue_gtp_annotated$variant_chrom <- factor(tissue_gtp_annotated$variant_chrom, levels = 1:22)

results_to_plot <- tissue_gtp_annotated[!is.na(chrom) & is_neand == TRUE]
setorder(results_to_plot, adj_p)
results_to_plot[, symbol_dedup := ""]
results_to_plot[!duplicated(symbol), symbol_dedup := symbol]

genes_to_highlight <- unique(results_to_plot[adj_p < 5e-8]$symbol)
genes_to_highlight <- genes_to_highlight[genes_to_highlight != ""]

don <- results_to_plot %>% 
  setorder(chrom, var_start) %>%
  
  # Compute chromosome size
  group_by(chrom) %>% 
  summarise(chr_len = max(var_start)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(results_to_plot, ., by = c("chrom" = "chrom")) %>%
  
  # Add a cumulative position of each SNP
  arrange(chrom, var_start) %>%
  mutate(BPcum = var_start + tot) %>%
  mutate(ital_symbol = paste0("italic(", symbol_dedup, ")"))

axisdf <- don %>% group_by(chrom) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

manhattan_plot <- ggplot(don, aes(x = BPcum, y = -log10(adj_p), label = symbol_dedup)) +
  
  # Show all points
  geom_point(aes(color = as.factor(chrom)), alpha = 0.8, size = 0.5) +
  scale_color_manual(values = rep(c("black", "grey"), 22 )) +
  geom_point(data = subset(don, symbol %in% genes_to_highlight), color = "red", size = 0.6) +
  
  # custom X axis:
  scale_x_continuous( label = axisdf$chrom, breaks = axisdf$center ) +
  scale_y_continuous(expand = c(0, 1)) +     # remove space between plot area and x axis
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid = element_blank()
  ) +
  geom_text_repel(data = subset(filter(don, adj_p < 5e-8), symbol_dedup != ""), fontface = "italic", size = 3.5, hjust = -0.1, vjust = 0.1, point.padding = NA) +
  xlab("Chromosome") +
  ylab(expression(-log[10](italic(p)))) +
  geom_hline(yintercept = -log10(5e-8), lty = "dashed", color = "gray")

plot_grid(plot_grid(ascertainment_plot, NULL, sqtl_count_plot, rel_widths = c(1.05, 0.1, 1), labels = c("A", "", "B"), nrow = 1), 
          NULL, 
          plot_grid(NULL, manhattan_plot, rel_widths = c(0.0175, 1)), 
          nrow = 3, align = "h", labels = c("", "", "C"), rel_heights = c(1, 0.05, 1), axis = "tb")

quartz.save("~/data/aseyedi2/manhattan/fig1.pdf", type = "pdf", height = 8, width = 8)