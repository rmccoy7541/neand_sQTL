#!/bin/bash
#SBATCH --job-name=sra2bam
#SBATCH --time=12:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-447%100

ml samtools
ml sra-tools

line=`sed "${SLURM_ARRAY_TASK_ID}q;d" sralist.txt`

sam-dump $line | samtools view -bS - > ${line}.bam
