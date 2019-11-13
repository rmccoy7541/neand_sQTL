library(tidyr)
library(dplyr)
library(data.table)

introns <- fread("../intron_counts/GTEx_v8_junctions_nohead.gct.gz",stringsAsFactors=FALSE,header=TRUE)
vcf <- fread("../vcf/vcf_for_merge.txt",stringsAsFactors=FALSE)
tissue_names <-c("Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Adrenal_Gland", "Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Brain_Amygdala", "Brain_Anterior_cingulate_cortex_BA24", "Brain_Caudate_basal_ganglia", "Brain_Cerebellar_Hemisphere", "Brain_Cerebellum", "Brain_Cortex", "Brain_Frontal_Cortex_BA9", "Brain_Hippocampus", "Brain_Hypothalamus", "Brain_Nucleus_accumbens_basal_ganglia", "Brain_Putamen_basal_ganglia", "Brain_Spinal_cord_cervical_c-1", "Brain_Substantia_nigra", "Breast_Mammary_Tissue", "Cells_Cultured_fibroblasts", "Cells_EBV-transformed_lymphocytes", "Colon_Sigmoid", "Colon_Transverse", "Esophagus_Gastroesophageal_Junction", "Esophagus_Mucosa", "Esophagus_Muscularis", "Heart_Atrial_Appendage", "Heart_Left_Ventricle", "Kidney_Cortex", "Liver", "Lung", "Minor_Salivary_Gland", "Muscle_Skeletal", "Nerve_Tibial", "Ovary", "Pancreas", "Pituitary", "Prostate", "Skin_Not_Sun_Exposed_Suprapubic", "Skin_Sun_Exposed_Lower_leg", "Small_Intestine_Terminal_Ileum", "Spleen", "Stomach", "Testis", "Thyroid", "Uterus", "Vagina", "Whole_Blood")

#get the paths for all sqtl perm pass results per tissue
tissue_names <- paste0("~/work/aseyedi2/sQTLv8/data/GTEx_Analysis_v8_sQTL/", tissue_names, ".v8.sqtl_signifpairs.txt.gz")

#read them all in
sqtl_by_tissue <- do.call(rbind, lapply(tissue_names, function(x) sqtl <- fread(x)))

sqtl_sep <- separate(sqtl_by_tissue, phenotype_id, c("chrom","start","end","cluster_id","ENSEMBL_ID"), sep=":", remove=TRUE)

sqtl_by_tissue <- NULL

combined <- inner_join(sqtl_sep, introns, by=c("ENSEMBL_ID"="Description"))

write.table(combined, file="combined.txt", sep="\t", quote=FALSE)

# join intron-sqtl and vcf dataframes by variant ID
# if there are multiple intron clusters that correspond to one variant,
# duplicate the line from the vcf
combined2 <- inner_join(combined, vcf, by=c("variant_id"="V1"))

write.table(combined2, file="combined2.txt", sep="\t", quote=FALSE)

# combine identical sample columns from vcf and intron files
# want to end up with a single sample column with info: GT;intron read count