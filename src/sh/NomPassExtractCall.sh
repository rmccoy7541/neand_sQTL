#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=NomPassExtractCall
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-100

######################
# Begin work section #
######################

ml R

line=`sed "${SLURM_ARRAY_TASK_ID}q;d" WHLBLD_chunks.txt`

Rscript --vanilla NomPassExtract.R ${line} "tag_snps.neand.EUR.bed"



