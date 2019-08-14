#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=NomPass
#SBATCH --nodes=1
#SBATCH --time=4:0:0
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --array=1-100

######################
# Begin work section #
######################

# nominal pass QTLtools
QTLtools cis \
  --vcf  $VCF \
  --bed "$pheno" \
  --cov  "$covariates" \
  --nominal 1  \
  --chunk $SLURM_ARRAY_TASK_ID 100 \
  --out "${tissue}_nominals_chunk_${SLURM_ARRAY_TASK_ID}.txt"

# extracts nom pass
sbatch -a $SLURM_ARRAY_TASK_ID \
  --export=listPath=$PWD/$abb,tissue=$(echo $abb),scripts=$scripts,worDir=$worDir,sprime=$sprime \
  ${scripts}sh/10_NomPassExtractCall.sh

# Call permuatation pass
sbatch -a $SLURM_ARRAY_TASK_ID \
  --export=VCF=$VCF,pheno=$(echo $abb/$abb.pheno.bed.gz),tissue=$(echo $abb/$abb),covariates=$(echo $abb/$full.v7.covariates_output.txt),abb=$abb,full=$full \
  ${scripts}/sh/11_PermPass.sh