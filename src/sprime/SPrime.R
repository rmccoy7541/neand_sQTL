library(data.table)
library(tidyverse)
library(parallel)

# this script reads sprime output and annotates matching (homozygous or heterozygous)
# to archaic genome VCFs. The tag "notcomp" refers to cases where a confident archaic genotype
# could not be assigned due to coverage, mapping, or other quality issues.

# TO DO: replace hard-coded paths with command-line arguments
args = commandArgs(trailingOnly=TRUE)
# args[1] is results.score, args[2] is splice ai directory
sprime <- fread(args[1])

compare_gt <- function(ref_allele, alt_allele, archaic_allele, genotype) {
  if (startsWith(genotype, "./.")) {
    return("notcomp") 
  } else if (startsWith(genotype, "0/0")) {
    if (ref_allele == archaic_allele) { 
      return("match") 
    } else {
      return("mismatch")
    }
  } else if (startsWith(genotype, "0/1")) {
    if (ref_allele == archaic_allele | alt_allele == archaic_allele) {
      return("match") 
    } else {
      return("mismatch")
    }
  } else if (startsWith(genotype, "1/1")) {
    if (alt_allele == archaic_allele) {
      return("match") 
    } else {
      return("mismatch")
    }
  } else {
    return("notcomp")
  }
}

get_archaic_gt <- function(sprime_line, archaic_dir) {
  
  # extract putative archaic introgressed allele from sprime output using ALLELE index
  introgressed_allele <- unlist(strsplit(paste(sprime_line$REF, sprime_line$ALT, sep = ","), split = ","))[sprime_line$ALLELE + 1]
  
  # get archaic genotypes from VCF
  archaic_vcf <- paste0(archaic_dir, "/chr", sprime_line$CHROM, "_mq25_mapab100.vcf.gz")
  archaic_gt <- suppressWarnings(fread(cmd = paste0("tabix ", archaic_vcf, " ", sprime_line$CHROM, ":", sprime_line$POS, "-", sprime_line$POS), sep = "\t"))
  
  # test whether archaic genotypes match putative introgressed allele
  if (nrow(archaic_gt) == 0) {
    return(data.table(sprime_line, altai_match = "notcomp", vindija_match = "notcomp", denisova_match = "notcomp"))
  } else {
    archaic_gt <- archaic_gt[, c(1:5, 10:12)] %>%
      setnames(., c("chrom", "pos", "id", "ref", "alt", "altai", "vindija", "denisova"))
   
    sprime_line[, altai_match := compare_gt(archaic_gt$ref, archaic_gt$alt, introgressed_allele, archaic_gt$altai)]
    sprime_line[, vindija_match := compare_gt(archaic_gt$ref, archaic_gt$alt, introgressed_allele, archaic_gt$vindija)]
    sprime_line[, denisova_match := compare_gt(archaic_gt$ref, archaic_gt$alt, introgressed_allele, archaic_gt$denisova)]
    
    return(sprime_line)
  }
}

# apply to all SNPs in parallel
results <- do.call(rbind, 
                   mclapply(1:nrow(sprime), 
                            function(x) get_archaic_gt(sprime[x,], args[2]), 
                            mc.cores = getOption("mc.cores", 8L)))

fwrite(results, 
       file = "sprime_calls.txt", 
       quote = FALSE, 
       sep = "\t", 
       col.names = TRUE, 
       row.names = FALSE)