## enrichment test on QTLs
library(data.table)
library(R.utils)
library(vcfR)
library(Matching)

cmdArgs = commandArgs(trailingOnly=TRUE)
# enrichment function
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


setwd(cmdArgs[1])
# gonna have to insert these things as cmd-line variables.
#setwd("Documents/GitHub/neand_sQTL/scratch/")

# 2 is vcf 
GTEx <- read.vcfR(cmdArgs[2], verbose = TRUE)
#GTEx <- read.vcfR("GTExSample.vcf", verbose = TRUE)
nomQTL <- fread(cmdArgs[3])
#nomQTL <- fread("THYROID_nominals_chunk_1_sample.txt")

QT_head <- c("phen.id", "phen.chrom.id", "phen.start", "phen.end", "strand.orient", "num.var.cis", "distance", "var.id",
             "var.chrom.id", "var.start", "var.end", "p.value", "reg.slope", "top.var")

nomQTL <- setNames(nomQTL, QT_head)

GTEx <- vcfR2tidy(GTEx, info_fields = "AF", info_types = c(AF = "n"))

GTEx <- data.table(GTEx[["fix"]]$CHROM, GTEx[["fix"]]$POS, GTEx[["fix"]]$ID, GTEx[["fix"]]$AF)

names(GTEx) <- c("var.chrom.id", "var.start", "var.id", "AF")

nomQTL <- subset(nomQTL, select = var.chrom.id:p.value)

GTEx <- transform(GTEx, var.chrom.id = as.integer(var.chrom.id))
AllNomFreq <- merge(nomQTL, GTEx, by = c("var.chrom.id", "var.start"))

permPass <- fread(cmdArgs[4])
#permPass <- fread("../analysis/PermPassResults/TopGenes_PermPass.txt")

# permHead <- c("phen.id", "phen.chrom.id", "phen.start", "phen.end", "strand.orient", "num.var.cis", "distance", "var.id",
#               "var.chrom.id", "var.start", "var.end", "deg.of.freedom.p.value", "dummy", "first.beta.dist", "second.beta.dist", 
#               "nom.p.value", "reg.slope", "empirical.p.value", "adjusted.p.value")
permPass <- setNames(permPass, c("intron.cluster", "var.id", "tissue.id", "adj.p", "q.val", "gene.sym"))

permPass <- merge(permPass, GTEx, by = "var.id")

enrichTest <- function(NLpermTest, AF) {
   
  
}