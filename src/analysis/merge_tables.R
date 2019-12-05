### USAGE ###
# Rscript table_merge_SY.R <intron counts file> <sQTL file> <tissue_id>
library(tidyr)
library(dplyr)
library(data.table)

# command line input specifies which intron counts and sQTL file to read
# and also the tissue name
args = commandArgs(trailingOnly=TRUE)

# read in intron counts file for one specific tissue
introns <- fread(args[2],stringsAsFactors=FALSE,header=TRUE)

# read in sQTL file for one specific tissue
sqtl <- fread(args[3], stringsAsFactors=FALSE, header=TRUE)

# separate intron cluster field to get ENSEMBL ID
sqtl_sep <- separate(sqtl, phenotype_id, c("chrom","start","end","cluster_id","ENSEMBL_ID"), sep=":", remove=TRUE)

# read in VCF file
vcf <- fread("../vcf/vcf_for_merge.txt.gz", stringsAsFactors=FALSE, header=TRUE)

tissue_name <- args[1]

dt <- inner_join(inner_join(sqtl_sep, introns, by=c("ENSEMBL_ID"="Description")), vcf, by=c("variant_id"="ID"))

### Get NL Isoforms

colnames(dt)[17] <- "transcript_id"

list_dt <- split.default(dt, nchar(names(dt)) > 10)

xcrips <- cbind(list_dt[[1]][,"variant_id"], list_dt[[2]][,"transcript_id"], list_dt[[2]][,5:ncol(list_dt[[2]])])

nl_iso <- cbind(list_dt[[1]][,"variant_id"], list_dt[[2]][,"transcript_id"], list_dt[[1]][,13:ncol(list_dt[[1]])])

#add variant id
xcrips <- as.data.table(xcrips %>% pivot_longer(-c(transcript_id, variant_id), names_to = "tissue_id", values_to = "counts"))

#add variant id
nl_iso <- as.data.table(nl_iso %>% pivot_longer(-c(transcript_id, variant_id), names_to = "individual", values_to = "is_NL"))

xcrips$individual <- gsub("^([^.]*.[^.]*)..*$", "\\1", xcrips$tissue_id)

final <- as.data.table(dplyr::full_join(xcrips, nl_iso, by = c("transcript_id", "individual", "variant_id")))

#final <- na.omit(final[is_NL == 1])

write.table(final,
            file = paste0(tissue_name, "_NL_isos.txt"),
            sep = "\t",
            row.names = F,
            quote = FALSE)