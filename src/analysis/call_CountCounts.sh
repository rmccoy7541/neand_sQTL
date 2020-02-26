#!/bin/bash
#SBATCH --job-name=CountCounts
#SBATCH --time=4:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=12
# number of tasks (processes) per node
#SBATCH --array=1-49
#######################################
# This is to call the script that counts the number of transcript/variant ID pairs
# with 0 counts for each genotype and with more than 0 counts for each genotype

ml R/3.5.1
ml gcc

line=`sed "${SLURM_ARRAY_TASK_ID}q;d" fileList.txt`

Rscript countCounts.R $line $line
