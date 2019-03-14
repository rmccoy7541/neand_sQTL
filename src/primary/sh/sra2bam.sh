#!/bin/bash
#SBATCH --job-name=sra2bam
#SBATCH --time=16:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=24
## YOU MUST ENTER ARRAY RANGE
#SBATCH --array=1-266%100

line=`sed "${SLURM_ARRAY_TASK_ID}q;d" ${sraListPath}/sralist.txt`

time {
	${homeDir}/progs/sra-tools/bin/sam-dump $line | ${homeDir}/progs/samtools-1.9/samtools view --threads 23 -bS - > ${line}.bam
}


