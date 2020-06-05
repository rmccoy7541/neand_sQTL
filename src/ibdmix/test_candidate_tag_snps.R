library(data.table)
library(pbmcapply)
library(tidyverse)
library(qqman)
library(rtracklayer)
library(liftOver)
library(Biostrings)

get_lod_stats <- function(chromosome, basedir, lod_threshold, length_threshold) {
  message(paste0("Reading data for chromosome ", chromosome, "."))
  
  file_path <- paste0(basedir, "gtex_ibdmix_chr", chromosome, ".txt.gz")
  chr_haps <- fread(file_path) %>%
    .[, length := end - start] %>%
    .[slod >= lod_threshold & length >= length_threshold]
  
  lod_stats <- data.table(pos = as.integer(unlist(strsplit(chr_haps$SNPs, ","))),
                          lod = as.numeric(unlist(strsplit(chr_haps$SNP_LODs, ",")))) %>%
    group_by(., pos) %>%
    summarize(., chrom = chromosome, mean_lod = mean(lod), max_lod = max(lod), sum_lod = sum(lod), n_haps = n()) %>%
    as.data.table()
  
  return(lod_stats)
}

lod_stat_summary <- do.call(rbind, pbmclapply(1:22, 
                                              function(x) get_lod_stats(x, "/work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/ibdmix_out/", 4, 50000), 
                                              mc.cores = getOption("mc.cores", 24L)))

tag_snp_query <- data.table(lod_stat_summary[, c("chrom", "pos")]) %>%
  .[, pos_1 := pos + 1]

invisible(lapply(1:22, function(x) fwrite(tag_snp_query[chrom == x], 
                                          file = paste0("/work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/tag_snps/tag_snp_query/tag_snp_query_chr", x, ".txt"),
                                          col.names = FALSE, row.names = FALSE, 
                                          sep = "\t", 
                                          quote = FALSE)))

retrieve_tag_snp_gt <- function(chromosome, wrkdir) {
  sys_command <- paste0("tabix ", wrkdir, "tag_snps/altai_gtex_combined_gt_chr", chromosome, ".bed.gz ",
                        "-R ", wrkdir, "tag_snps/tag_snp_query/tag_snp_query_chr", chromosome, ".txt ", 
                        "> ", wrkdir, "tag_snps/tag_snp_gt_chr", chromosome, ".txt")
  system(sys_command)
}

invisible(pbmclapply(1:22, 
                     function(x) retrieve_tag_snp_gt(x, "/work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/"), 
                     mc.cores = getOption("mc.cores", 24L)))

get_gt <- function(chromosome, basedir) {
  message(paste0("Reading genotype data for chromosome ", chromosome, "."))
  
  file_path <- paste0(basedir, "tag_snp_gt_chr", chromosome, ".txt")
  chr_gt <- fread(file_path) %>%
    .[, -ncol(.), with = FALSE]
  
  header <- fread(paste0(basedir, "header.txt"), header = FALSE)
  setnames(chr_gt, header$V1)
  
  return(chr_gt)
}

get_introgressed_haplotypes <- function(chromosome, basedir, lod_threshold = 4, length_threshold = 50000) {
  message(paste0("Reading data for chromosome ", chromosome, "."))
  
  file_path <- paste0(basedir, "gtex_ibdmix_chr", chromosome, ".txt.gz")
  chr_haps <- fread(file_path) %>%
    .[, length := end - start] %>%
    .[slod >= lod_threshold & length >= length_threshold]
  
  return(chr_haps)
}

introgressed_haplotypes <- do.call(rbind, pbmclapply(1:22, 
                                                     function(x) get_introgressed_haplotypes(x, "/work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/ibdmix_out/"), 
                                                     mc.cores = getOption("mc.cores", 24L)))

fisher_test_snp <- function(genotypes, haplotypes, lod_stats, index) {
  
  query_chrom <- lod_stats[index,]$chrom
  query_pos <- lod_stats[index,]$pos
  
  hap_subset <- haplotypes[chrom == query_chrom & start <= query_pos & end > query_pos]
  introgressed_samples <- hap_subset$ID
  
  genotypes_snp <- genotypes[chrom == query_chrom & pos == query_pos][1,]
  
  gt_subset <- genotypes_snp %>%
    .[, c(6:ncol(.)), with = FALSE] %>%
    t() %>%
    as.data.table(keep.rownames = TRUE) %>%
    setnames(., c("sample", "gt")) %>%
    .[gt != 9]
  
  ref_allele <- genotypes_snp$ref
  alt_allele <- genotypes_snp$alt
  
  ac <- sum(gt_subset[-1,]$gt)
  an <- length(gt_subset[-1,]$gt) * 2
  af <- ac / an
  
  introgressed_an <- nrow(gt_subset[sample %in% introgressed_samples]) * 2
  if (length(introgressed_samples) == 0) {
    introgressed_ac <- 0
  } else {
    introgressed_ac <- sum(gt_subset[sample %in% introgressed_samples]$gt)
  }
  nonintrogressed_an <- nrow(gt_subset[!(sample %in% introgressed_samples)]) * 2
  nonintrogressed_ac <- sum(gt_subset[!(sample %in% introgressed_samples)]$gt)
  
  fisher_results <- fisher.test(rbind(cbind(introgressed_ac, 
                                            introgressed_an - introgressed_ac), 
                                      cbind(nonintrogressed_ac, 
                                            nonintrogressed_an - nonintrogressed_ac)))
  
  return(data.table(ref = ref_allele, alt = alt_allele, af = af,
                    int_ac = introgressed_ac, int_an = introgressed_an,
                    non_ac = nonintrogressed_ac, non_an = nonintrogressed_an,
                    fisher_p = fisher_results$p.value))
    
}

pbmclapply_with_error <- function(X,FUN,...){    
  pbmclapply(X, function(x, ...) tryCatch(FUN(x, ...),
                                      error=function(e) NULL),
             mc.cores = getOption("mc.cores", 24L))
}

fisher_test_wrapper <- function(chromosome, basedir, haplotypes, lod_stats, cores) {
  
  message(paste0("Computing Fisher's exact test for candidate archaic SNPs on chromosome ", chromosome, "."))
  
  lod_stats_subset <- lod_stats[chrom == chromosome]
  haplotypes_subset <- haplotypes[chrom == chromosome]
  
  chr_genotypes <- get_gt(chromosome, basedir)
  
  fisher <- rbindlist(pbmclapply_with_error(1:nrow(lod_stats_subset), 
                                            function(x) fisher_test_snp(chr_genotypes, haplotypes_subset, lod_stats_subset, x), 
                                            mc.cores = getOption("mc.cores", cores)))
  
  lod_stats_subset <- cbind(lod_stats_subset, fisher)
  
  chr_genotypes[, chr_pos := paste(chrom, pos, sep = "_")]
  lod_stats_subset[, chr_pos := paste(chrom, pos, sep = "_")]
  lod_stats_subset <- merge(lod_stats_subset, chr_genotypes[, c("chr_pos", "AltaiNea_AltaiNea")], "chr_pos")
  
  return(lod_stats_subset)
  
}

get_pan_tro <- function(chain, fasta_file, input_snps, index) {
  input_snp <- input_snps[index,] %>%
    .[, chr := paste0("chr", chrom)] %>%
    .[, strand := "+"] %>%
    makeGRangesFromDataFrame(seqnames.field = "chr", start.field = "pos", end.field = "pos")
    
  seqlevelsStyle(input_snp) <- "UCSC"
  
  chimp_coords <- liftOver(input_snp, chain) %>%
    as.data.table() %>%
    .[1,]
  
  if (nrow(chimp_coords) == 0) {
    return(NA)
  }
  
  chimp_query <- paste0(chimp_coords$seqnames, ":", chimp_coords$start, "-", chimp_coords$end)
  chimp_strand <- chimp_coords$strand
  
  chimp_allele <- toupper(fread(cmd = paste("/home-net/home-4/rmccoy22@jhu.edu/code/samtools-1.10/bin/samtools faidx", 
                                            fasta_file, chimp_query, sep = " "))[[1]])
  
  if (chimp_strand == "-") {
    rev_comp <- reverseComplement(DNAString(chimp_allele)) %>%
      as.character()
    return(rev_comp)
  } else if (chimp_strand == "+") {
    return(chimp_allele)
  }
}

tag_snp_stats <- rbindlist(lapply(1:22, 
                                  function(x) fisher_test_wrapper(x, "/work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/tag_snps/tag_snp_gt/", 
                                                                  introgressed_haplotypes, lod_stat_summary, 24L)))


candidate_archaic_snps <- tag_snp_stats[fisher_p < 5e-8]

chain_input <- import.chain("/work-zfs/rmccoy22/resources/panTro6/hg38.panTro6.all.chain")

fasta_path <- "/work-zfs/rmccoy22/resources/panTro6/panTro6.fa.gz"

chimp_alleles <- pbmclapply(1:nrow(candidate_archaic_snps),
                            function(x) tryCatch({get_pan_tro(chain = chain_input, 
                                                              fasta_file = fasta_path, 
                                                              input_snps = candidate_archaic_snps, 
                                                              index = x)}, 
                                                 error = function(e) NA), 
                            mc.cores = getOption("mc.cores", 24L))

candidate_archaic_snps[, chimp_allele := unlist(chimp_alleles)]
candidate_archaic_snps[, neand_1 := as.character(NA)]
candidate_archaic_snps[, neand_2 := as.character(NA)]
candidate_archaic_snps[AltaiNea_AltaiNea == 0, neand_1 := ref]
candidate_archaic_snps[AltaiNea_AltaiNea == 0, neand_2 := ref]
candidate_archaic_snps[AltaiNea_AltaiNea == 1, neand_1 := ref]
candidate_archaic_snps[AltaiNea_AltaiNea == 1, neand_2 := alt]
candidate_archaic_snps[AltaiNea_AltaiNea == 2, neand_1 := alt]
candidate_archaic_snps[AltaiNea_AltaiNea == 2, neand_2 := alt]


candidate_archaic_snps[, neand_is_derived := !is.na(chimp_allele) &
                         (chimp_allele == ref | chimp_allele == alt) & 
                         (neand_1 != chimp_allele | neand_2 != chimp_allele)]

candidate_archaic_snps[, is_archaic_snp := FALSE]
candidate_archaic_snps[neand_1 == alt & neand_2 == alt &
                         neand_is_derived == TRUE & fisher_p < 5e-8 & mean_lod > 0 &
                         (int_ac / int_an) > (non_ac / non_an),
                       is_archaic_snp := TRUE]
candidate_archaic_snps[neand_1 == ref & neand_2 == ref &
                         neand_is_derived == TRUE & fisher_p < 5e-8 & mean_lod > 0 &
                         (int_ac / int_an) < (non_ac / non_an),
                       is_archaic_snp := TRUE]

ggplot(data = candidate_archaic_snps[is_archaic_snp == TRUE], aes(x = daf)) +
  geom_histogram() +
  theme_classic()
