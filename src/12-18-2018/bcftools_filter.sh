#!/bin/bash
#SBATCH --job-name=bialleliconly
#SBATCH --time=8:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=24

ml bcftools
# filters all non-biallelic sites
bcftools view -m2 -M2 -v snps --threads 23 -O z -o GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz ${1}/ncbi/files/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz