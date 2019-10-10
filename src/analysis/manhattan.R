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
library(rtracklayer)

count_sqtl <- function(tissue, summarize = FALSE, select_neand = TRUE) {
  # gtp <- fread(paste0(snakemake@input[[1]]))
  # gtp <- fread("/scratch/groups/rmccoy22/aseyedi2/sQTLv8/data/GTEx_Analysis_v8_sQTL/Whole_Blood_permutation_table_NE.txt")
  gtp <- fread(paste0(tissue, "_permutation_table_NE.txt")) 
  
  # neand <- fread(snakemake@input[[2]])[vindija_match == "match" | altai_match == "match"] %>%
  neand <- fread("~/data/aseyedi2/neanderthal-sqtl/metadata/sprime_calls.txt")[vindija_match == "match" | altai_match == "match"] %>%
    mutate(., var_id_1 = paste(CHROM, POS, REF, ALT, "b38", sep = "_")) %>%
    as.data.table()
  
  gtp[, is_neand := variant_id %in% neand$var_id]
  
  gtp_neand_t <- gtp[is_neand == TRUE]
  gtp_neand_f <- gtp[is_neand == FALSE]
  
  gtp_neand_t[, logP := -log10(pval_nominal)]
  setorder(gtp_neand_t, logP)
  gtp_neand_t[, expectedP := rev(-log10(ppoints(n = length(gtp_neand_t$pval_nominal))))]
  
  gtp_neand_f[, logP := -log10(pval_nominal)]
  setorder(gtp_neand_f, logP)
  gtp_neand_f[, expectedP := rev(-log10(ppoints(n = length(gtp_neand_f$pval_nominal))))]
  
  gtp <- rbind(gtp_neand_t, gtp_neand_f)
  
  if (summarize == TRUE) {
    return(data.table(TISSUE = tissue, n_sqtl = nrow(gtp[is_neand == select_neand]), n_total = nrow(gtp)))
  } else {
    if (count(gtp) > 0) {
      print(tissue)
      return(data.table(gtp, TISSUE_ID = tissue))
    }
    else {
      print(paste0("Tissue ", tissue, " has no samples with NL introgression."))
      return(NULL)
    }
  }
}

tissue_names <-c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra", "Breast_Mammary_Tissue", "Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Kidney_Cortex", "Liver", "Lung", "Minor_Salivary_Gland", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")
tissue_abbv <- c("ADPSBQ", "ADPVC", "ADRNLG", "ARTAORT", "ARTCRN", "ARTTBL", "BRNAMY", "BRNACC", "BRNCDT", "BRNCHB", "BRNCHA", "BRNCTXA", "BRNCTXB", "BRNHPP", "BRNHPT", "BRNNCC", "BRNPTM", "BRNSPC", "BRNSNG", "BREAST", "FIBRBLS", "LCL", "CLNSGM", "CNLTRN", "ESPGEJ", "ESPMCS", "ESPMSL", "HRTAA", "HRTLV", "KDNCTX", "LIVER", "LUNG", "SLVRYG", "MSCLSK", "NERVET", "OVARY", "PNCREAS", "PTTARY", "PRSTTE", "SKINNS", "SKINS", "SNTTRM", "SPLEEN", "STMACH", "TESTIS", "THYROID", "UTERUS", "VAGINA", "WHLBLD")

names(tissue_names) <- tissue_abbv

sqtl_counts_by_tissue <- do.call(rbind, lapply(tissue_names, function(x) count_sqtl(x, summarize=TRUE)))
sqtl_counts_by_tissue[, n_samples := c(763,564,275,450,253,770,177,213,291,263,298,325,425,243,236,277,232,182,164,480,527,192,389,432,401,622,559,452,689,100,251,867,181,1132,722,195,360,301,262,638,849,193,260,381,406,812,166,173,3288)]

tissue_table <- data.table(TISSUE_ID = tissue_names, SAMPLE_SIZE = c(763,564,275,450,253,770,177,213,291,263,298,325,425,243,236,277,232,182,164,480,527,192,389,432,401,622,559,452,689,100,251,867,181,1132,722,195,360,301,262,638,849,193,260,381,406,812,166,173,3288))

tissue_gtp <- do.call(rbind, lapply(tissue_names, function(x) count_sqtl(x)))
tissue_gtp <- merge(tissue_gtp, tissue_table, "TISSUE_ID")

#gene_list <- rtracklayer::import(snakemake@input["gtf"]) %>%
#/scratch/groups/rmccoy22/aseyedi2/sQTLv8/data/GTEx_Analysis_v8_sQTL

dplyr::select(head(tissue_gtp[is_neand == TRUE], 100), phenotype_id, variant_id, TISSUE_ID, pval_nominal, gene_name)

# tissue_gtp$variant_chrom <- factor(tissue_gtp$variant_chrom, levels = 1:22)

results_to_plot <- tissue_gtp[is_neand == TRUE]
setorder(results_to_plot, pval_nominal)
results_to_plot[, symbol_dedup := ""]
results_to_plot[!duplicated(gene_name), symbol_dedup := gene_name]
results_to_plot$seqnames <- as.integer(gsub('[a-zA-Z]', '', results_to_plot$seqnames))
colnames(results_to_plot)[16] <- "chrom"

genes_to_highlight <- unique(results_to_plot[pval_nominal < 5e-8]$gene_name)
genes_to_highlight <- genes_to_highlight[genes_to_highlight != ""]


# TODO: Figure out the deal with var_start. There is a column named "start" but that's hg37. 


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