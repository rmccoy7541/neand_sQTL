#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=PermPass
#SBATCH --nodes=1
#SBATCH --time=16:0:0
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3

######################
# Begin work section #
######################

/home-1/aseyedi2@jhu.edu/work/progs/QTLtools/QTLtools_1.1_Ubuntu14.04_x86_64 cis \
  --vcf  $VCF \
  --bed "$pheno" \
  --cov  "${covariates}" \
  --permute 1000  \
  --chunk ${SLURM_ARRAY_TASK_ID} 100 \
  --out "${tissue}_permutations_chunk_${SLURM_ARRAY_TASK_ID}.txt" \
  --normal

scripts=$(echo /home-1/aseyedi2@jhu.edu/work/aseyedi2/neand_sQTL/src/primary/)

sbatch -a $SLURM_ARRAY_TASK_ID --export=VCF=$VCF,pheno=$(echo $abb/$abb.pheno.bed.gz),tissue=$(echo $abb/$abb),covariates=$(echo $abb/$full.v7.covariates_output.txt),permutations=$(echo $abb/$abb.permutations_full_FDR.thresholds.txt) ${scripts}/sh/CondPass.sh

sbatch --export=abb=$(echo $tissue) $scripts/../AWS/PermPassFDR.sh
