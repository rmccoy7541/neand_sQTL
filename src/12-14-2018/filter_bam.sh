#!/bin/bash
#SBATCH --job-name=filter_bam
#SBATCH --time=3:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-10

# load samtools
ml samtools

# bamlist.txt is list of bam files generated from sra's using `ls >> bamlist.txt`
line=`sed "${SLURM_ARRAY_TASK_ID}q;d" bamlist.txt`

filt="${line}.filt"

echo "Filtering $line"

samtools view -L GRCh37.bed $line > "$filt"
