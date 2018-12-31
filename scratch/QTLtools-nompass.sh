#!/bin/bash
#SBATCH --job-name=QTLtools-nompass
#SBATCH --time=2:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-20

./QTLtools_1.1_Ubuntu14.04_x86_64 cis --vcf GTExWholeGenomeSequenceGenotypeMatrixBiallelicOnly.vcf.gz --bed testNE_sQTL_perind.counts.gz.qqnorm_chr${SLURM_ARRAY_TASK_ID}.gz --cov Whole_Blood.v7.covariates_output.txt --nominal 0.01 --chunk 1 1 --out nominals_${SLURM_ARRAY_TASK_ID}.txt
