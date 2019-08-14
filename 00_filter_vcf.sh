#!/usr/bin/env bash
# This script assumes that the genotype file has already been downloaded. It filters non-biallelic sites and indexes it.

## Step 4 - Genotype & Covariates Preparation
################################################
# insert scripts dir below
scripts=$(echo {path/to/project/scripts})
# insert the full path for GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz
vcf=$(echo {/path/to/}phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz)
outdir=$(echo {desired/output/directory})

# filter the non-biallelic sites from genotype file using bcftools; see script for details
sbatch --wait --export=vcf=$vcf,outdir=$outdir ${scripts}/sh/00a_bcftools_filter.sh
# index our friend with tabix
echo "Indexing our friend..."
sbatch --export=outdir=$outdir $scripts/sh/00b_index_vcf.sh