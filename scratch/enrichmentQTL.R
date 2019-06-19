## enrichment test on QTLs
library(data.table)
library(R.utils)
# enrichment funciton
# - We are finding random sets of non-NL introgressed SNPs that have a similar frequency to the NL-introgressed SNPs and 
# determining how many of them are significant to build a null distribution. We then take the signficant NL-introgressed SNPs 
# and determine where they lie along the distribution. Permuations are the testing of the random set of SNPs.
# 
# From nominal pass column, take `uniq` - that SNP will have an associated frequency
# 
# From the VCF, grab all the SNPs that are in the nominal pass file so then you'll have frequency annotated with NL from SPrime
# 
# do this for every permutation to find frequency matched SNPs
# 
# Use `pymatch` or `Matching` for frequency matching
# 
# Develop seed method for random permutation - set default seed

setwd("Documents/GitHub/neand_sQTL/scratch/")

vcf <- fread("GTExVCFSample.vcf.gz")
nomQTL <- fread("THYROID_nominals_chunk_1_sample.txt")

QT_head <- c("phen.id", "phen.chrom.id", "phen.start", "phen.end", "strand.orient", "num.var.cis", "distance", "var.id",
             "var.chrom.id", "var.start", "var.end", "p-value", "reg.slope", "top.var")

nomQTL <- setNames(nomQTL, QT_head)




