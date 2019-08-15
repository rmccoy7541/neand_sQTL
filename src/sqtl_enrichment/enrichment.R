library(data.table)
library(tidyverse)
library(broom)
library(Matching)
library(pbmcapply)
library(qvalue)

# cmd_args[1] is tissue name, e.g., THYROID
cmd_args <- commandArgs(trailingOnly = TRUE)
tissue_input <- cmd_args[1]
# tissue_input <- "THYROID"

# read perm pass
perm <- fread(paste0("/work-zfs/rmccoy22/aseyedi2/sqtl_permutation_backup/", tissue_input, "_permutations.txt")) %>%
  setnames(., c("intron_cluster", "chrom", "pheno_start", "pheno_end", 
                "strand", "total_cis", "distance", "variant_id", "variant_chrom", 
                "var_start", "var_end", "df", "dummy", "param_1", "param_2",
                "p", "beta", "emp_p", "adj_p"))

perm[, qval := qvalue(perm$adj_p)$qvalues]
perm[, is_sig := qval < 0.1]

# read Neanderthal SNP annotations
sprime <- fread("/work-zfs/rmccoy22/aseyedi2/neanderthal-sqtl/analysis/SPRIME/sprime_calls.txt")
sprime[, variant_id := paste(CHROM, POS, REF, ALT, "b37", sep = "_")]
sprime[, is_neand := (altai_match == "match" | vindija_match == "match")]
conf_neand_snps <- sprime[is_neand == TRUE & !(grepl(",", variant_id))]$variant_id
lowc_neand_snps <- sprime[is_neand == FALSE & !(grepl(",", variant_id))]$variant_id

# remove Sprime-significant, but non-matching SNPs, as their introgression status is uncertain
perm <- perm[!(variant_id %in% lowc_neand_snps)]
perm[, is_neand := variant_id %in% conf_neand_snps]

# get tissue-specific allele frequencies
command <- paste("grep -v ^#", paste0("/work-zfs/rmccoy22/rmccoy22/sqtl/enrichment/tissue_vcf/", tissue_input, ".vcf"), "| cut -f3,7-8 | tr ';' '\t' | cut -f1-2,4 | sed 's/AF=//g'")
af <- fread(cmd = command, header = FALSE) %>%
  setnames(., c("variant_id", "FILTER", "AF"))
# limit to stringent QC "PASS" SNPs
af <- af[FILTER == "PASS"]
af[, AF := as.numeric(AF)]
af[, MAF := AF]
af[AF > 0.5, MAF := 1 - AF]

perm <- merge(perm, af, "variant_id")

# get ld scores
ldsc <- fread("/work-zfs/rmccoy22/rmccoy22/sqtl/enrichment/ldsc/l2.ldscore.gz", header = TRUE) %>%
  setnames(., c("CHR", "variant_id", "BP", "L2"))
ldsc <- ldsc[, c("variant_id", "L2"), with = FALSE]

perm <- merge(perm, ldsc, "variant_id")

rm(ldsc)
rm(af)

# naive enrichment test
fisher_no_covar <- tidy(fisher.test(
  rbind(
    cbind(nrow(perm[is_neand == FALSE & is_sig == FALSE]), nrow(perm[is_neand == FALSE & is_sig == TRUE])),
    cbind(nrow(perm[is_neand == TRUE & is_sig == FALSE]), nrow(perm[is_neand == TRUE & is_sig == TRUE]))
  )
)) %>%
  as.data.table()
fisher_no_covar[, method := "no_covar"]


# matched enrichment test
# function: shuffles permutation-pass sQTL data, matches Neand / non-Neand SNPs on covariates
matched_enrichment <- function(input_data) {
  shuffled_input <- input_data[sample(nrow(input_data), replace = FALSE),]
  # match on MAF and LD score
  matches <- Match(X = as.matrix(shuffled_input[, c("MAF", "L2")]), Tr = shuffled_input$is_neand, caliper = 0.1, replace = FALSE, ties = FALSE, version = "fast", M = 100)
  output_data <- rbind(shuffled_input[unique(matches$index.control),], shuffled_input[unique(matches$index.treated),])
  fisher_results <- fisher.test(
    rbind(
      cbind(nrow(output_data[is_neand == FALSE & is_sig == FALSE]), nrow(output_data[is_neand == FALSE & is_sig == TRUE])),
      cbind(nrow(output_data[is_neand == TRUE & is_sig == FALSE]), nrow(output_data[is_neand == TRUE & is_sig == TRUE]))
    )
  )
  return(fisher_results)
}

RNGkind("L'Ecuyer-CMRG")
set.seed(123)
mc.reset.stream()
fisher_covar <- tidy(matched_enrichment(perm)) %>%
  as.data.table()
fisher_covar[, method := "covar"]

fwrite(rbind(fisher_no_covar, fisher_covar), 
  paste0("/work-zfs/rmccoy22/rmccoy22/sqtl/enrichment/output/", tissue_input, ".txt"), 
  quote = FALSE, row.names = FALSE, col.names = FALSE)
