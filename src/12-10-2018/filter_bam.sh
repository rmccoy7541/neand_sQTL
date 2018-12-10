#!/bin/bash
#SBATCH --job-name=filter_bam
#SBATCH --time=1:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-10

# load samtools
ml samtools

# create text file with all *.bam files
# ls *.bam | tail -n +2 >> bamlist.txt
# run this as a one-off command; job array will create massive list with duplicates

# each line with a .SRA file is a job in the job array
line=`sed "${SLURM_ARRAY_TASK_ID}q;d" bamlist.txt`

filt="${line}.filt"

echo "Filtering $line"

samtools view -b -L GRCh37.bed $line > "$filt"
