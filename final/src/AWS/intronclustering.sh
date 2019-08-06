#!/bin/bash
#SBATCH --job-name=intronclustering
#SBATCH --time=6:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=12

ml python/2.7-anaconda

mkdir intronclustering/

python /scratch/groups/rmccoy22/aseyedi2/leafcutter/clustering/leafcutter_cluster.py -j juncfiles.txt -r intronclustering/ -m 50 -o Ne-sQTL -l 500000