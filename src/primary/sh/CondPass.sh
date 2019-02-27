#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=CondPass
#SBATCH --nodes=1
#SBATCH --time=5:0:0
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-100

######################
# Begin work section #
######################

./QTLtools_1.1_Ubuntu14.04_x86_64 cis \
	--vcf  $VCF \
	--bed "$pheno" \
	--cov  "Whole_Blood.v7.covariates_output.txt" \
	--mapping "permuatations_full_FDR.thresholds.txt" \
	--chunk ${SLURM_ARRAY_TASK_ID} 100 \
	--out "WHLBLD_permutations_chunk_${SLURM_ARRAY_TASK_ID}.txt"
