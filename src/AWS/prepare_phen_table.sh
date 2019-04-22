#!/bin/bash
#SBATCH --job-name=prepare_phen_table
#SBATCH --time=2:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=6


ml python/2.7-anaconda
python /scratch/groups/rmccoy22/aseyedi2/leafcutter/scripts/prepare_phenotype_table.py Ne-sQTL_perind.counts.gz -p 10