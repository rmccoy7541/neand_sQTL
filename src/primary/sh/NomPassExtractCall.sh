#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=NomPassExtractCall
#SBATCH --time=4:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-100%50

######################
# Begin work section #
######################

module load R
module load gcc

# Consider passing WHLBLD_chunks as a command line var
line=`sed "${SLURM_ARRAY_TASK_ID}q;d" ${listPath}/${tissue}_chunks.txt`

Rscript --vanilla ../R/NomPassExtract.R ${line} "tag_snps.neand.EUR.bed"



