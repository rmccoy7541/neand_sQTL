library(data.table)
library(tidyverse)

af <- fread("gtex_af.txt") %>%
  setnames(., c("variant_id", "freq"))

nom <- fread("THYROID_nominals_allchunks.txt") %>%
  setnames(., c("phen.id", "phen.chrom.id", "phen.start", "phen.end", "strand.orient", "num.var.cis", "distance", "var.id",
                "var.chrom.id", "var.start", "var.end", "p.value", "reg.slope", "top.var"))

# modify below to read in permutation results
perm <- fread("Permutationpass") %>%
  setnames(., c("phen.id", "phen.chrom.id", "phen.start", "phen.end", "strand.orient", "num.var.cis", "distance", "var.id",
                "var.chrom.id", "var.start", "var.end", "deg.of.freedom.p.value", "dummy", "first.beta.dist", "second.beta.dist", 
                "nom.p.value", "reg.slope", "empirical.p.value", "adjusted.p.value"))

# edit below to read in sprime Neanderthal SNPs 
sprime_snps <- fread("sprime_file.txt")$variant_id

# restrict frequency data to SNPs that were actually tested for associations
af <- af[variant_id %in% nom$var.id]

# annotate each snp as Neanderthal or not using same definition as we used previously
af[, is_neand := variant_id %in% sprime_snps]

# annotate each snp as signficant or not based on permutation results
af[, is_sig := variant_id %in% perm[q < 0.01]$var.id]

run_permutation <- function(< list of inputs here >) {
  # code here to sample random SNPs with (is_neand == FALSE)
  # and count how many are significant based on permutation pass
  return(< results of one permutation >)
}

null_distribution <- do.call(rbind, lapply(1:1000, function(x) run_permutation())
                             
 # then compare real data to distribution of significant hits under null