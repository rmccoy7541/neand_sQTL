#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=CondPass
#SBATCH --nodes=1
#SBATCH --time=1:0:0
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --array=1-100

######################
# Begin work section #
######################

/home-1/aseyedi2@jhu.edu/work/progs/QTLtools/QTLtools_1.1_Ubuntu14.04_x86_64 cis \
	--vcf  $VCF \
	--bed "$pheno" \
	--cov  "$covariates" \
	--mapping "$permutations" \
	--chunk ${SLURM_ARRAY_TASK_ID} 100 \
	--out "${tissue}_conditionals_chunk_${SLURM_ARRAY_TASK_ID}.txt"
