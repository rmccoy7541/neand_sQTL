library(data.table)
library(tidyverse)
library(qvalue)
library(pbmcapply)
library(Matching)

# cmd_args[1] is tissue name - "THYROID"
# cmd_args[2] is allele freq - "/work-zfs/rmccoy22/aseyedi2/GTExWGS_VCF/GTExWGS.AF.all.txt"
# cmd_args[3] is base path - "/work-zfs/rmccoy22/aseyedi2/sqtl_permutation_backup/all_noms/varIDs/chunks/"
# cmd_args[4] is sprime - "/work-zfs/rmccoy22/aseyedi2/neanderthal-sqtl/analysis/SPRIME/sprime_calls.txt"
# cmd_args[5] is perm pass file - "/work-zfs/rmccoy22/aseyedi2/sqtl_permutation_backup/THYROID_permutations.txt"
# cmd_args[6] is output file 
# cmd_args[7] is seed (def 123)
# cmd_args[8] is M (1000?)

##cmd_args <- commandArgs(trailingOnly = TRUE)
tissue_input <- cmd_args[1]
#tissue_input <- "THYROID"

# allele freq VCF
#af <- fread(file = "/work-zfs/rmccoy22/aseyedi2/GTExWGS_VCF/GTExWGS.AF.all.txt", colClasses = c("character", "integer", "character", "character", "numeric"))
af <- fread(file = cmd_args[2], colClasses = c("character", "integer", "character", "character", "numeric"))
af[, variant_id := paste(CHROM, POS, REF, ALT, "b37", sep = "_")]
af[, c("CHROM", "POS", "REF", "ALT") := NULL]

# compute minor allele frequency
af[, MAF := AF]
af[AF > 0.5, MAF := 1 - AF]

basedir <- cmd_args[3]
#basedir <- "/work-zfs/rmccoy22/aseyedi2/sqtl_permutation_backup/all_noms/varIDs/chunks/"

noms <- unique(fread(paste0("/work-zfs/rmccoy22/rmccoy22/sqtl/nom_varID_chunks_uniq/", 
                            tissue_input, "_nom_varIDs_uniq.txt"), header = FALSE)$V1)

# subset allele frequencies to those tested for splicing associations (nominal pass)
af <- af[variant_id %in% noms]

# NL SNPs
sprime <- fread(cmd_args[4])
#sprime <- fread("/work-zfs/rmccoy22/aseyedi2/neanderthal-sqtl/analysis/SPRIME/sprime_calls.txt")

sprime[, variant_id := paste(CHROM, POS, REF, ALT, "b37", sep = "_")]

# create column for if either match
sprime[, is_neand := (altai_match == "match" | vindija_match == "match")]

neand_snps <- sprime[is_neand == TRUE & !(grepl(",", variant_id))]$variant_id

af[, is_neand := variant_id %in% neand_snps]

# read perm pass
#perm <- fread(paste0("/work-zfs/rmccoy22/aseyedi2/sqtl_permutation_backup/", tissue_input, "_permutations.txt")) %>%
perm <- fread(cmd_args[5]) %>%
  setnames(., c("intron_cluster", "chrom", "pheno_start", "pheno_end", 
                "strand", "total_cis", "distance", "variant_id", "variant_chrom", 
                "var_start", "var_end", "df", "dummy", "param_1", "param_2",
                "p", "beta", "emp_p", "adj_p"))

# finding q value; FDR
perm[, qval := qvalue(perm$adj_p)$qvalues]
sig_snps <- perm[qval < 0.1]$variant_id

# record significance
af[, is_sig := variant_id %in% sig_snps]

setorder(af, MAF)

# create an allele frequency bin factor variable
af[, MAF_bin := factor(round(MAF, 3))]

# compute partial contigency tables, stratifying on allele frequency bin
partial_tables <- xtabs(~ is_sig + is_neand + MAF_bin, data = af)

# convert the tables to numeric rather than integer (to prevent overflow)
partial_tables_numeric <- data.table(apply(partial_tables, 2, as.numeric))
partial_tables_numeric[, MAF_bin := rep(as.numeric(unique(af$MAF_bin)), each = 2)]

partial_tables_array <- simplify2array(by(partial_tables_numeric, partial_tables_numeric$MAF_bin, as.matrix))[, -3,]

mh_results <- mantelhaen.test(partial_tables_array)

output_data <- data.table(
  tissue = tissue_input, 
  est = mh_results$estimate,
  ci.ll = mh_results$conf.int[1],
  ci.ul = mh_results$conf.int[2],
  p = mh_results$p.value
)

#fwrite(output_data, file = "test.out", quote = F, col.names = T, row.names = F)
fwrite(output_data, file = cmd_args[6], quote = F, col.names = T, row.names = F)