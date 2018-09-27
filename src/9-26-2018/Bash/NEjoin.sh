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

join -1 1 -2 1 tag_snps.neand.EUR.sorted.txt rsIDchr$SLURM_ARRAY_TASK_ID.txt | tr ' ' '\t' >> tag_snps.neand.EUR.rsid.txt
join -1 1 -2 1 tag_snps.neand.EAS.sorted.txt rsIDchr$SLURM_ARRAY_TASK_ID.txt | tr ' ' '\t' >> tag_snps.neand.EAS.rsid.txt
join -1 1 -2 1 tag_snps.neand.SAS.sorted.txt rsIDchr$SLURM_ARRAY_TASK_ID.txt | tr ' ' '\t' >> tag_snps.neand.SAS.rsid.txt