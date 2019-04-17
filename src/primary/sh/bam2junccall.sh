#!/bin/bash
#SBATCH --job-name=bam2junccall
#SBATCH --time=4:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=1

# creates list of filtered bam files for job array
# ls *.filt >> filtbamlist.txt

# each line is a filtered bam file to be used as an input in the job array
line=`sed "${SLURM_ARRAY_TASK_ID}q;d" filtlist.txt`

# load samtools
# ml samtools
ml python/2.7-anaconda

# converts each file to junc
# Remember to set LeafCutter full path (leafCutterDir) in bam2junc.sh

# IMPORTANT: changed bam2junc.sh in leafcutter/scripts/ to directly call the
# bin file for samtools 
echo "Converting ${line} to junc"
sh /scratch/groups/rmccoy22/aseyedi2/leafcutter/scripts/bam2junc.sh ${line} ${line}.junc

