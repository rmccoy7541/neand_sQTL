#!/bin/bash
#SBATCH --job-name=bam2junccall
#SBATCH --time=1:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-10

# creates list of filtered bam files for job array
ls *.filt | tail -n +2 >> filtbamlist.txt

# each line is a bam file to be used as an input in the job array
line=`sed "${SLURM_ARRAY_TASK_ID}q;d" filtbamlist.txt`

# load samtools
ml samtools

# converts each file to junc
# Remember to set LeafCutter full path (leafCutterDir) in bam2junc.sh
echo "Converting ${line} to ${line}.junc"
sh ../scripts/bam2junc.sh ${line} ${line}.junc
echo ${line}.junc >> test_juncfiles.txt
