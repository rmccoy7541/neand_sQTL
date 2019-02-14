#!/bin/bash
#SBATCH --job-name=QTLtools-NomPass
#SBATCH --time=4:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-20

# tissue is going to be $tissue

./QTLtools_1.1_Ubuntu14.04_x86_64 cis \
  --vcf  GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz \
  --bed "$tissue" \
  --cov  "Whole_Blood.v7.covariates_output.txt" \
  --nominal 0.01 \
  --chunk $SLURM_ARRAY_TASK_ID 20 \
  --out "$WHLBLD_nominals_chunk_${SLURM_ARRAY_TASK_ID}.txt"
