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

cmd_args <- commandArgs(trailingOnly = TRUE)
tissue_input <- cmd_args[1]
# tissue_input <- "THYROID"

# allele freq VCF
# af <- fread(file = "/work-zfs/rmccoy22/aseyedi2/GTExWGS_VCF/GTExWGS.AF.all.txt", colClasses = c("character", "integer", "character", "character", "numeric"))
af <- fread(file = cmd_args[2], colClasses = c("character", "integer", "character", "character", "numeric"))
af[, AF := as.numeric(AF)]
af[, variant_id := paste(CHROM, POS, REF, ALT, "b37", sep = "_")]
af[, c("CHROM", "POS", "REF", "ALT") := NULL]

basedir <- cmd_args[3]
# basedir <- "/work-zfs/rmccoy22/aseyedi2/sqtl_permutation_backup/all_noms/varIDs/chunks/"

# collects all chunks of computed nominal pass output
read_nom_ids_wrapper <- function(basedir, tissue_name) {
  basepath <- paste0(basedir, tissue_name, "_nom_varIDs_chunk")
  dt <- do.call(c, pbmclapply(sprintf("%03d", 1:100), function(x) read_nom_ids(basepath, x), mc.cores = getOption("mc.cores", 24L)))
  return(dt)
}

# reads all nominal pass chunks
read_nom_ids <- function(basepath, chunk_number) {
  path <- paste0(basepath, "_", chunk_number, ".txt")
  return(unique(fread(path, select = 1, header = F)$V1))
}

# calls the above two functions
noms <- unique(read_nom_ids_wrapper(basedir, tissue_input))

# subset allele frequencies to those tested for splicing associations (nominal pass)
af <- af[variant_id %in% noms]

# NL SNPs
sprime <- fread(cmd_args[4])
# sprime <- fread("/work-zfs/rmccoy22/aseyedi2/neanderthal-sqtl/analysis/SPRIME/sprime_calls.txt")

sprime[, variant_id := paste(CHROM, POS, REF, ALT, "b37", sep = "_")]

# create column for if either match
sprime[, is_neand := (altai_match == "match" | vindija_match == "match")]

neand_snps <- sprime[is_neand == TRUE & !(grepl(",", variant_id))]$variant_id

af[, is_neand := variant_id %in% neand_snps]

# read perm pass
# perm <- fread("/work-zfs/rmccoy22/aseyedi2/sqtl_permutation_backup/THYROID_permutations.txt") %>%
perm <- fread(cmd_args[5]) %>%
  setnames(., c("intron_cluster", "chrom", "pheno_start", "pheno_end", 
                "strand", "total_cis", "distance", "variant_id", "variant_chrom", 
                "var_start", "var_end", "df", "dummy", "param_1", "param_2",
                "p", "beta", "emp_p", "adj_p"))

# finding q value; FDR
perm[, qval := qvalue(perm$adj_p)$qvalues]

sig_snps <- perm[qval < 0.1]$variant_id

af[, is_sig := variant_id %in% sig_snps]

# shuffle AF
# set.seed(123)
set.seed(cmd_args[7])
af <- af[sample(1:nrow(af), replace = F),]

# play around with changing "M"?
matches <- Match(X = af$AF, Tr = af$is_neand, exact = TRUE, replace = FALSE, ties = FALSE, version = "fast", M = cmd_args[8])

# enrichment test

sig_neand <- nrow(af[unique(matches$index.treated),][is_sig == TRUE,])
nsig_neand <- nrow(af[unique(matches$index.treated),][is_sig == FALSE,])
sig_nneand <- nrow(af[unique(matches$index.control),][is_sig == TRUE,])
nsig_nneand <- nrow(af[unique(matches$index.control),][is_sig == FALSE,])

print(paste0("Significant NL: ", sig_neand))
print(paste0("Insignificant NL: ", nsig_neand))
print(paste0("Significant non-NL:", sig_nneand))
print(paste0("Non-significant non-NL: ", nsig_nneand))

fisher_results <- fisher.test(rbind(cbind(sig_neand, nsig_neand), cbind(sig_nneand, nsig_nneand)))

output_data <- data.table(
  tissue = tissue_input, 
  est = fisher_results$estimate,
  ci.ll = fisher_results$conf.int[1],
  ci.ul = fisher_results$conf.int[2],
  p = fisher_results$p.value
  )
# fwrite(output_data, file = "test.out", quote = F, col.names = T, row.names = F)
fwrite(output_data, file = cmd_args[6], quote = F, col.names = T, row.names = F)