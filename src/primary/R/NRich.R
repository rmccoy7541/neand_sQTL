library(data.table)
library(tidyverse)
library(qvalue)
library(Matching)

af <- fread(file = "~/Downloads/GTExWGS.AF.txt", 
            colClasses = c("character", "integer", "character", "character", "numeric"))

af[, AF := as.numeric(AF)]

sprime <- fread("~/Desktop/sQTL/sprime_calls.txt.gz")

af[, variant_id := paste(CHROM, POS, REF, ALT, "b37", sep = "_")]
sprime[, variant_id := paste(CHROM, POS, REF, ALT, "b37", sep = "_")]

af[, c("CHROM", "POS", "REF", "ALT") := NULL]

sprime[, is_neand := (altai_match == "match" | vindija_match == "match")]

neand_snps <- sprime[is_neand == TRUE & !(grepl(",", variant_id))]$variant_id

af[, is_neand := variant_id %in% neand_snps]

perm <- fread("~/Desktop/sQTL/ADPSBQ_permutations.txt") %>%
  setnames(., c("intron_cluster", "chrom", "pheno_start", "pheno_end", 
                "strand", "total_cis", "distance", "variant_id", "variant_chrom", 
                "var_start", "var_end", "df", "dummy", "param_1", "param_2",
                "p", "beta", "emp_p", "adj_p"))

perm[, qval := qvalue(perm$adj_p)$qvalues]

sig_snps <- perm[qval < 0.1]$variant_id

af[, is_sig := variant_id %in% sig_snps]

# naive enrichment test

fisher.test(rbind(table(af[is_neand == FALSE]$is_sig),
                  table(af[is_neand == TRUE]$is_sig))
           )