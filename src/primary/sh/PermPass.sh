#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=PermPass
#SBATCH --nodes=1
#SBATCH --time=3:0:0
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --array=1-100

######################
# Begin work section #
######################

/home-1/aseyedi2@jhu.edu/work/progs/QTLtools/QTLtools_1.1_Ubuntu14.04_x86_64 cis \
  --vcf  GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz \
  --bed "WHLBLD.pheno.bed.gz" \
  --cov  "Whole_Blood.v7.covariates_output.txt" \
  --permute 1000  \
  --chunk ${SLURM_ARRAY_TASK_ID} 100 \
  --out "WHLBLD_permutations_chunk_${SLURM_ARRAY_TASK_ID}.txt"
  --normal
