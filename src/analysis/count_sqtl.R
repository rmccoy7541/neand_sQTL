library(data.table)
library(tidyverse)
library(pbapply)
library(ggrepel)
library(annotate)
library(rtracklayer)

count_sqtl <- function(tissue, summarize = FALSE) {
  gtp <- fread(paste0(tissue, "_permutation_table_NE.txt")) 

  gtp[, logP := -log10(pval_nominal)]
  setorder(gtp, logP)
  gtp[, expectedP := rev(-log10(ppoints(n = length(gtp$pval_nominal))))]

  if (summarize == TRUE) {
    return(data.table(TISSUE = tissue, n_sqtl = nrow(gtp)))
  } 
  else if (nrow(gtp) < 2) {
    gtp <- NULL
    tissue <- NULL
    }
  else {
    return(data.table(gtp, TISSUE_ID = tissue))
  }
}

tissue_names <-c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra", "Breast_Mammary_Tissue", "Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Kidney_Cortex", "Liver", "Lung", "Minor_Salivary_Gland", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")
tissue_abbv <- c("ADPSBQ", "ADPVC", "ADRNLG", "ARTAORT", "ARTCRN", "ARTTBL", "BRNAMY", "BRNACC", "BRNCDT", "BRNCHB", "BRNCHA", "BRNCTXA", "BRNCTXB", "BRNHPP", "BRNHPT", "BRNNCC", "BRNPTM", "BRNSPC", "BRNSNG", "BREAST", "FIBRBLS", "LCL", "CLNSGM", "CNLTRN", "ESPGEJ", "ESPMCS", "ESPMSL", "HRTAA", "HRTLV", "KDNCTX", "LIVER", "LUNG", "SLVRYG", "MSCLSK", "NERVET", "OVARY", "PNCREAS", "PTTARY", "PRSTTE", "SKINNS", "SKINS", "SNTTRM", "SPLEEN", "STMACH", "TESTIS", "THYROID", "UTERUS", "VAGINA", "WHLBLD")

names(tissue_names) <- tissue_abbv

sqtl_counts_by_tissue <- do.call(rbind, lapply(tissue_names, function(x) count_sqtl(x, summarize=TRUE)))
sqtl_counts_by_tissue[, n_samples := c(581,469,233,387,213,584,129,147,194,175,209,205,175,165,170,202,170,126,114,396,483,147,318,368,330,497,465,372,386,73,208,515,144,706,532,167,305,237,221,517,605,174,227,324,322,574,129,141,670)]




### number of sQTLs per tissue
# png(filename = "sQTLs_per_tissue.png")
ggplot(data = sqtl_counts_by_tissue, aes(x = n_samples, y = n_sqtl, label = names(tissue_names), color = TISSUE)) +
  theme_bw() +
  geom_point() +
  ggrepel::geom_text_repel(force = 1) + # hjust = 0, nudge_x = 4, nudge_y = 2
  coord_cartesian(xlim = c(0, 900)) +
  coord_cartesian(ylim = c(0,115)) +
  theme(panel.grid = element_blank(), legend.position = "none") +
  xlab("Number of Genotyped Samples") +
  ylab("Number of sQTLs")
# dev.off()
###

tissue_gtp <-  do.call(rbind, lapply(tissue_names, function(x) count_sqtl(x, summarize = FALSE))) %>%
  setorder(., pval_nominal)

topGenes <- dplyr::select(tissue_gtp, phenotype_id, variant_id, TISSUE_ID, pval_nominal, pval_beta, gene_name) %>%
  unique(.)

write.csv(topGenes, file = "TopGenes_PermPass_All.csv", eol = "\r\n", row.names = F)
