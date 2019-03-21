#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=NomPass
#SBATCH --nodes=1
#SBATCH --time=12:0:0
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --array=1-100

######################
# Begin work section #
######################

/scratch/groups/rmccoy22/progs/QTLtools/QTLtools_1.1_Ubuntu14.04_x86_64 cis \
	--vcf  $VCF \
	--bed "$pheno" \
	--cov  "$covariates" \
	--nominal 1  \
	--chunk $SLURM_ARRAY_TASK_ID 100 \
	--out "${tissue}_nominals_chunk_${SLURM_ARRAY_TASK_ID}.txt"
