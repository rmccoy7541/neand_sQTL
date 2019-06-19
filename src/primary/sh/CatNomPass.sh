#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=CatNominals
#SBATCH --nodes=1
#SBATCH --time=30:0:0
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --array=1-48
######################
# Begin work section #
######################

cd /work-zfs/rmccoy22/rmccoy22/sqtl/intron_clustering/tryagain

TISSUE=`sed "${SLURM_ARRAY_TASK_ID}q;d" tissuesused.txt`

for i in {1..100}; do
    echo "Catting ${TISSUE}_nominals_chunk_${i}.txt..."
    cat ${TISSUE}/${TISSUE}_nominals_chunk_${i}.txt >> ${TISSUE}_nominals.txt
done

mv ${TISSUE}_nominals.txt $TISSUE