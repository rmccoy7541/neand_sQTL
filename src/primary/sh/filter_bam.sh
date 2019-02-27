#!/bin/bash
#SBATCH --job-name=filter_bam
#SBATCH --time=4:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=24
#SBATCH --array=1-354%50

# bamlist.txt is list of bam files generated from sra's using `ls >> bamlist.txt`
line=`sed "${SLURM_ARRAY_TASK_ID}q;d" bamlist.txt`

filt="${line}.filt"

echo "Filtering $line"

/scratch/groups/rmccoy22/progs/samtools-1.9/samtools view --threads 23 -b -L GRCh37.bed $line > "$filt"


