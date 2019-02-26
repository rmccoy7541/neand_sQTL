#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=NomPassExtractCall
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-100%50

######################
# Begin work section #
######################

ml R

# Consider passing WHLBLD_chunks as a command line var
line=`sed "${SLURM_ARRAY_TASK_ID}q;d" ${listPath}/WHLBLD_chunks.txt`

Rscript NomPassExtract.R ${line} "tag_snps.neand.EUR.bed"



