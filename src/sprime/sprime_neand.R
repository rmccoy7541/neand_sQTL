library(data.table)
library(tidyverse)
library(parallel)
library(pbmcapply)
library(rtracklayer)

# this script reads sprime output and annotates matching (homozygous or heterozygous)
# to archaic genome VCFs. The tag "notcomp" refers to cases where a confident archaic genotype
# could not be assigned due to coverage, mapping, or other quality issues.

args = commandArgs(trailingOnly = TRUE)
# snakemake@input[[1]] is results.score file (the output of SPrime), 
# snakemake@input[[2]] is directory containing merged archaic VCFs ("/scratch/users/rmccoy22@jhu.edu/archaic_splicing/spliceai/")
sprime <- fread(snakemake@input[[1]])

# liftover the coordinates to hg19 to retrieve archaic alleles from their corresponding VCFs
sprime[, strand_orientation := "+"]
hg38_coords <- makeGRangesFromDataFrame(sprime, 
                                        seqnames.field = "CHROM", 
                                        start.field = "POS", 
                                        end.field = "POS", 
                                        strand.field = "strand_orientation", 
                                        ignore.strand = FALSE)

chain <- import.chain("/work-zfs/rmccoy22/resources/reference/liftover/hg38ToHg19.over.chain")
hg19_coords <- liftOver(hg38_coords, chain) %>%
  as.data.table() %>%
  setnames(., "group", "index")
sprime[, index := .I]
sprime <- merge(sprime, hg19_coords, "index", all.x = TRUE)

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
  archaic_vcf <- paste0(archaic_dir, "/", sprime_line$seqnames, "_mq25_mapab100.vcf.gz")
  if (is.na(sprime_line$seqnames)) {
    return(data.table(sprime_line, altai_match = "notcomp", vindija_match = "notcomp", denisova_match = "notcomp"))
  }
  archaic_gt <- suppressWarnings(fread(cmd = paste0("tabix ", archaic_vcf, " ", gsub("chr", "", sprime_line$seqnames), ":", sprime_line$start, "-", sprime_line$start), sep = "\t"))
  
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
                   pbmclapply(1:nrow(sprime), 
                              function(x) get_archaic_gt(sprime[x,], snakemake@input[[2]]), 
                              mc.cores = getOption("mc.cores", 48L)))

fwrite(results, 
       file = "sprime_calls.txt", 
       quote = FALSE, 
       sep = "\t", 
       col.names = TRUE, 
       row.names = FALSE)