#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=PermPass
#SBATCH --nodes=1
#SBATCH --time=10:0:0
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=5

######################
# Begin work section #
######################

/home-1/aseyedi2@jhu.edu/work/progs/QTLtools/QTLtools_1.1_Ubuntu14.04_x86_64 cis \
  --vcf  $VCF \
  --bed $tissue.pheno.bed.gz \
  --cov  "$covariates" \
  --permute 1000  \
  --chunk $SLURM_ARRAY_TASK_ID 100 \
  --out "${tissue}_permutations_chunk_${SLURM_ARRAY_TASK_ID}.txt" \
  --normal
  
 #sbatch -a 1-100 --export=VCF=$VCF,pheno=$(echo $abb/$abb.pheno.bed.gz),tissue=$(echo $abb/$abb),covariates=$(echo $abb/$full.v7.covariates_output.txt),abb=$abb,full=$full ${scripts}/sh/PermPass.sh