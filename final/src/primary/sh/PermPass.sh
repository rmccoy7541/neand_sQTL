#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=PermPass
#SBATCH --nodes=1
#SBATCH --time=30:0:0
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4

######################
# Begin work section #
######################

echo "Permutation Pass variables are: "
echo "Tissue - $tissue"
echo "Covariates - $covariates"
echo "Chunk - $SLURM_ARRAY_TASK_ID"

/home-1/aseyedi2@jhu.edu/work/progs/QTLtools/QTLtools_1.1_Ubuntu14.04_x86_64 cis \
  --vcf  $VCF \
  --bed $tissue.pheno.bed.gz \
  --cov  "$covariates" \
  --permute 1000  \
  --chunk $SLURM_ARRAY_TASK_ID 100 \
  --out "${tissue}_permutations_chunk_${SLURM_ARRAY_TASK_ID}.txt" \
  --normal
