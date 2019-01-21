#!/bin/bash
#SBATCH --job-name=QTLtools-NomPass
#SBATCH --time=2:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-22

# tissue is going to be $1 and chunk $2

./QTLtools_1.1_Ubuntu14.04_x86_64 cis \
  --vcf ../GTExWholeGenomeSequenceGenotypeMatrixBiallelicOnly.vcf.gz \
  --bed $tissue \
  --cov  \
  --nominal 0.01 \
  --chunk $chunk 20 \
  --out ${tissue}_nominals_chunk${chunk}.txt

