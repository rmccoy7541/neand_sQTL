#!/bin/bash
#SBATCH --job-name=sra2bam
#SBATCH --time=24:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=24
#SBATCH --array=1-454%100

ml samtools
ml sra-tools

line=`sed "${SLURM_ARRAY_TASK_ID}q;d" sralist.txt`

time {
	sam-dump $line | samtools view --threads 23 -bS - > ${line}.bam
}