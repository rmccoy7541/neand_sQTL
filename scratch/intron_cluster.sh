#!/bin/bash
#SBATCH --job-name=intron_cluster
#SBATCH --time=2:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-10

ml python/2.7-anaconda

python ../clustering/leafcutter_cluster.py -j ${SLURM_ARRAY_TASK_ID}.split -r intronclustering/ -m 50 -o testNE_sQTL -l 500000

# rm ${SLURM_ARRAY_TASK_ID}.split
