#!/bin/bash
#SBATCH --job-name=rename_gtex
#SBATCH --time=0:15:0
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-10

# load samtools
ml samtools

# bamlist.txt is list of bam files generated from sra's using `ls >> bamlist.txt`
file=`sed "${SLURM_ARRAY_TASK_ID}q;d" filtlist.txt`
# get GTEX ID from within file
ID=$(samtools view -H $file | grep -o -m 1 'GTEX-....')
# rename file
mv $file ${ID}
