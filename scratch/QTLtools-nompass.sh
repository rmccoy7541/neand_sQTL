#!/bin/bash
#SBATCH --job-name=QTLtools-nompass
#SBATCH --time=2:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-22

for tissue in testcovmerge/*; do ./QTLtools_1.1_Ubuntu14.04_x86_64 cis --vcf GTExWholeGenomeSequenceGenotypeMatrixBiallelicOnly.vcf.gz 
--bed 
testNE_sQTL_perind.counts.gz.qqnorm_chr${SLURM_ARRAY_TASK_ID}.gz --cov ${tissue} --nominal 0.01 --chunk 1 1 --out 
${tissue}_nominals_chr${SLURM_ARRAY_TASK_ID}.txt; done;