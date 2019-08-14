#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=NomPassExtractCall
#SBATCH --time=4:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4

######################
# Begin work section #
######################

module load R
module load gcc


line=$(echo ${tissue}_nominals_chunk_${SLURM_ARRAY_TASK_ID})

# extracts NL-introgressed nominal pass calls
Rscript ${scripts}/R/12_NomPassExtract.R $listPath/${line} $sprime