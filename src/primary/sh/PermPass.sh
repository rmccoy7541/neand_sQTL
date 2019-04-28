#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=PermPass
#SBATCH --nodes=1
#SBATCH --time=10:0:0
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4

######################
# Begin work section #
######################

/home-1/aseyedi2@jhu.edu/work/progs/QTLtools/QTLtools_1.1_Ubuntu14.04_x86_64 cis \
  --vcf  $VCF \
  --bed LUNG/LUNG.pheno.bed.gz \
  --cov  "LUNG/Lung.v7.covariates_output.txt" \
  --permute 1000  \
  --chunk 19 100 \
  --out "LUNG/LUNG_permutations_chunk_19.txt" \
  --normal
  
  sbatch -a 19 --export=VCF=$VCF,pheno=$(echo LUNG/LUNG.pheno.bed.gz),tissue=$(echo LUNG/LUNG),covariates=$(echo LUNG/Lung.v7.covariates_output.txt),abb=$(echo LUNG),full=$(echo Lung) ${scripts}/sh/PermPass.sh