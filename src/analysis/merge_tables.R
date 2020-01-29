### USAGE ###
# Rscript table_merge_SY.R <intron counts file> <sQTL file> <tissue_id>
library(tidyverse)
library(data.table)

# command line input specifies which intron counts and sQTL file to read
# and also the tissue name
args = commandArgs(trailingOnly=TRUE)

# read in intron counts file for one specific tissue
# introns <- fread(snakemake@input["introns"],stringsAsFactors=FALSE,header=TRUE)
# introns <- fread("Muscle_Skeletal_intronCounts.txt",stringsAsFactors=FALSE,header=TRUE)
# introns <- fread("../splitIC/Muscle_Skeletal_intronCounts.txt",stringsAsFactors=FALSE,header=TRUE)
introns <- fread(args[1],stringsAsFactors=FALSE,header=TRUE)


# read in sQTL file for one specific tissue
# sqtl <- fread(snakemake@input["perm"], stringsAsFactors=FALSE, header=TRUE)
# sqtl <- fread("/scratch/groups/rmccoy22/aseyedi2/sQTLv8/data/GTEx_Analysis_v8_sQTL/Muscle_Skeletal.v8.sqtl_signifpairs.txt.gz", stringsAsFactors=FALSE, header=TRUE)
# sqtl <- fread("Muscle_Skeletal.v8.sqtl_signifpairs.txt.gz", stringsAsFactors=FALSE, header=TRUE)
sqtl <- fread(args[2], stringsAsFactors=FALSE, header=TRUE)

# separate intron cluster field to get ENSEMBL ID
sqtl_sep <- separate(sqtl, phenotype_id, c("chrom","start","end","cluster_id","ENSEMBL_ID"), sep=":", remove=TRUE)

# read in VCF file
# vcf <- fread(snakemake@input["vcf_merge"], stringsAsFactors=FALSE, header=TRUE)
vcf <- fread("../vcf/vcf_for_merge.txt.gz", stringsAsFactors=FALSE, header=TRUE)

# tissue_name <- snakemake@input["tisname"]
tissue_name <- "Muscle_Skeletal"

dt <- inner_join(inner_join(sqtl_sep, introns, by=c("ENSEMBL_ID"="Description")), vcf, by=c("variant_id"="ID"))

### Get NL Isoforms

colnames(dt)[17] <- "transcript_id"

list_dt <- split.default(dt, nchar(names(dt)) > 10)

variant_id <- as.character(list_dt[[1]][,"variant_id"])
transcript_id <- as.character(list_dt[[2]][,"transcript_id"])

xcrips <- data.table(variant_id, list_dt[[2]][,5:ncol(list_dt[[2]])])

nl_iso <- data.table(variant_id, transcript_id, list_dt[[1]][,13:ncol(list_dt[[1]])])

rm(variant_id, transcript_id)

#add variant id
xcrips <- as.data.table(xcrips %>% pivot_longer(-c(transcript_id, variant_id), names_to = "tissue_id", values_to = "counts"))

#add variant id
nl_iso <- as.data.table(nl_iso %>% pivot_longer(-c(transcript_id, variant_id), names_to = "individual", values_to = "is_NL"))

xcrips$individual <- gsub("^([^-]*-[^-]*)-.*$", "\\1", xcrips$tissue_id)

# Joins tables, gets rid of all NAs and consolidates counts by matching xcrips, vars and individuals, and counts number of rows consolidated per match
final <- as.data.table(dplyr::full_join(xcrips, nl_iso, by = c("transcript_id", "variant_id", "individual"))) %>% 
  na.omit() %>%
  group_by(variant_id, transcript_id, is_NL) %>%
  dplyr::summarise(counts=sum(counts), nrows=n()) %>%
  as.data.table()

write.table(final,
            file = paste0(tissue_name, "_NL_isos.txt"),
            sep = "\t",
            row.names = F,
            quote = FALSE)
