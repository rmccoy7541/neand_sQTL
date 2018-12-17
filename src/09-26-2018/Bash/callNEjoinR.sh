#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=NEjoin
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-user=aseyedi2@jhu.edu
#SBATCH --array=1-22

######################
# Begin work section #
######################

module load gcc/5.5.0
module load R # default is 3.5
module list
which R # should give absolute path

Rscript --vanilla NEjoin.R tag_snps.neand.EAS.sorted.txt rsIDchr$SLURM_ARRAY_TASK_ID.txt $SLURM_ARRAY_TASK_ID 'EAS'
Rscript --vanilla NEjoin.R tag_snps.neand.EUR.sorted.txt rsIDchr$SLURM_ARRAY_TASK_ID.txt $SLURM_ARRAY_TASK_ID 'EUR'
Rscript --vanilla NEjoin.R tag_snps.neand.SAS.sorted.txt rsIDchr$SLURM_ARRAY_TASK_ID.txt $SLURM_ARRAY_TASK_ID 'SAS'
