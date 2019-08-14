#!/bin/bash
#SBATCH --job-name=QTLtools-Filter
#SBATCH --time=2:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --array=1-22


ml htslib

# get the phenotype files
line=`sed "${SLURM_ARRAY_TASK_ID}q;d" leafcutterphenotypes.txt`
# convert
cat ${line} | awk '{ $4=$4" . +"; print $0 }' | tr " " "\t" | bgzip -c > $line.qtltools
