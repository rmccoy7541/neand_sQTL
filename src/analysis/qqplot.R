library(data.table)
library(tidyverse)
library(ggrepel)
library(annotate)
library(org.Hs.eg.db)
library(rtracklayer)
library(pbapply)
library(parallel)

count_sqtl <- function(tissue, summarize = FALSE) {
    gtp <- fread(paste0(tissue, ".v8.sqtl_allpairs.txt.gz")) %>%
        setorder(., pval_nominal)
    
    gtp[, qval := qvalue(gtp$pval_nominal)$qvalues]
    
    neand <- fread("~/work/aseyedi2/neand_sqtl/metadata/sprime_calls.txt")[vindija_match == "match" | altai_match == "match"] %>%
        mutate(., var_id = paste(CHROM, POS, REF, ALT, "b37", sep = "_")) %>%
        as.data.table()
    
    gtp[, logP := -log10(pval_nominal)]
    setorder(gtp, logP)
    gtp[, expectedP := rev(-log10(ppoints(n = length(gtp$pval_nominal))))]
    
    gtp[, is_neand := variant_id %in% neand$var_id]
    
    if (summarize == TRUE) {
        return(data.table(TISSUE = tissue, n_sqtl = nrow(gtp[is_neand == TRUE & qval < 0.1]), n_total = nrow(gtp[qval < 0.1])))
    } else {
        print(paste("Reading ", tissue, "..."))
        return(data.table(gtp, TISSUE_ID = tissue))
    }
    
}

tissue_names <-c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra", "Breast_Mammary_Tissue", "Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Kidney_Cortex", "Liver", "Lung", "Minor_Salivary_Gland", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")
tissue_abbv <- c("ADPSBQ", "ADPVC", "ADRNLG", "ARTAORT", "ARTCRN", "ARTTBL", "BRNAMY", "BRNACC", "BRNCDT", "BRNCHB", "BRNCHA", "BRNCTXA", "BRNCTXB", "BRNHPP", "BRNHPT", "BRNNCC", "BRNPTM", "BRNSPC", "BRNSNG", "BREAST", "FIBRBLS", "LCL", "CLNSGM", "CNLTRN", "ESPGEJ", "ESPMCS", "ESPMSL", "HRTAA", "HRTLV", "KDNCTX", "LIVER", "LUNG", "SLVRYG", "MSCLSK", "NERVET", "OVARY", "PNCREAS", "PTTARY", "PRSTTE", "SKINNS", "SKINS", "SNTTRM", "SPLEEN", "STMACH", "TESTIS", "THYROID", "UTERUS", "VAGINA", "WHLBLD")

names(tissue_names) <- tissue_abbv

tissue_gtp <-  do.call(rbind, mclapply(tissue_names, function(x) count_sqtl(x, summarize = FALSE), mc.cores = 24L)) %>%
    setorder(., pval_nominal)

ggplot(data = tissue_gtp[is_neand == TRUE], aes(x = expectedP, y = logP, color = TISSUE_ID)) +
    theme_bw() +
    theme(panel.grid = element_blank(), legend.position = "none") +
    geom_point() +
    geom_abline(slope = 1, intercept = 0) +
    xlab(expression(Expected -log[10](italic("p")))) +
    ylab(expression(Observed -log[10](italic("p")))) +
    facet_wrap(~ TISSUE_ID)