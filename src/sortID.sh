#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=sortID
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=aseyedi2@jhu.edu
#SBATCH --array=1-22

######################
# Begin work section #
######################

awk '{print $1"_"$2"\t"$3}' rsID_chr$SLURM_ARRAY_TASK_ID.txt > rsIDchr$SLURM_ARRAY_TASK_ID.txt

