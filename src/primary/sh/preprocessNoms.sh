#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=preprocessNoms
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --array=1-48

######################
# Begin work section #
######################

line=`sed "${SLURM_ARRAY_TASK_ID}q;d" tissuenames.txt`

echo $line

for i in {1..100}; do
	cut -d' ' -f 8 /work-zfs/rmccoy22/rmccoy22/sqtl/intron_clustering/tryagain/$line/${line}_nominals_chunk_${i}.txt >> /work-zfs/rmccoy22/aseyedi2/sqtl_permutation_backup/all_noms/varIDs/${line}_nom_varIDs.txt
done
