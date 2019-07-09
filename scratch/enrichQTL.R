library(tidyverse)
library(data.table)
library(ff)

# 1 is the VCF filtered for allele freq, 2 is nom file, 3 is perm pass
cmdArgs = commandArgs(trailingOnly=TRUE)

# af <- fread(cmdArgs[1])
af <- fread("GTExWGS.AF.txt")

# nom <- fread("../sqtl_permutation_backup/all_noms/THYROID_nominals.txt")
nom <- fread("../sqtl_permutation_backup/all_noms/THYROID_nominals.txt") %>%
  setnames(., c("phen.id", "phen.chrom.id", "phen.start", "phen.end", "strand.orient", "num.var.cis", "distance", "var.id",
                "var.chrom.id", "var.start", "var.end", "p.value", "reg.slope", "top.var"))

# modify below to read in permutation results
perm <- fread("../sqtl_permutation_backup/THYROID_permutations.txt") %>%
  setnames(., c("phen.id", "phen.chrom.id", "phen.start", "phen.end", "strand.orient", "num.var.cis", "distance", "var.id",
                "var.chrom.id", "var.start", "var.end", "deg.of.freedom.p.value", "dummy", "first.beta.dist", "second.beta.dist", 
                "nom.p.value", "reg.slope", "empirical.p.value", "adjusted.p.value"))

# edit below to read in sprime Neanderthal SNPs 
sprime_snps <- fread("sprime_file.txt")$ID

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