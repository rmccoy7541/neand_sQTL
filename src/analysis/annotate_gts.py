#!/usr/bin/env python3

"""
Usage: ./annotate_gts.py

This script adds and annotates a field, ARCHAIC_GT, in the FORMAT field of each sample in a VCF. ARCHAIC_GT 
indicates whether the sample carries an archaic introgressed allele at that locus. The location of 
introgressed alleles and their archaic/modern identity were discovered with SPrime.

0 - homozygous modern human allele
1 - heterozygous for modern human/archaic allele
2 - homozygous for archaic allele
3 - unknown GT (i.e. sample was not genotyped in the VCF)
"""

import pysam

# sprime_allele = 0: ref allele is archaic
# sprime_allele = 1, 2, etc.: alt allele is archaic
def allele_match(GT, sprime_allele):
    matches = 0
    if GT[0] == sprime_allele:
        matches += 1
    if GT[1] == sprime_allele:
        matches += 1
    if GT[0] == None:
        matches = 3
    # 0 = homozygous MH, 1 = heterozygous MH/archaic, 2 = homozygous archaic, 3 = unknown GT
    return(matches)

# read in VCF file
vcf = pysam.VariantFile("GTEx_v8_SPrimeAnn.vcf.gz", "r")
# add ARCHAIC_GT field to each sample
vcf.header.formats.add("ARCHAIC_GT", ".", "Integer", "Match to archaic/modern allele")

# create VCF to write to
w = pysam.VariantFile("GTEx_v8_ArchaicGTAnn.vcf", "w", header=vcf.header)

# iterate through variants, annotating each sample with whether its genotype matches the archaic allele
with open("GTEx_v8_ArchaicGTAnn.vcf", "a") as outfile:
    for variant in vcf:
        # extract archaic allele from the variant's INFO field
        archaic_allele = variant.info["ARCHAIC_ALLELE"]
        for sample in variant.samples:
            # for each sample, return the sample's genotype for the archaic allele
            # 0 = homozygous MH, 1 = heterozygous MH/archaic, 2 = homozygous archaic, 3 = unknown GT
            archaic_gt = allele_match(variant.samples[sample]["GT"], archaic_allele)
            variant.samples[sample]["ARCHAIC_GT"] = archaic_gt
        
        outfile.write(str(variant))
