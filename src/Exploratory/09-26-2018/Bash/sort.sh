#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=sort
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=aseyedi2@jhu.edu
#SBATCH --array=1-22

######################
# Begin work section #
######################

sort -k 1,1 -o rsID_chr$SLURM_ARRAY_TASK_ID.txt rsID_chr$SLURM_ARRAY_TASK_ID.txt