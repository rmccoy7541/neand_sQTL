#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=NRich.sh
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=18
#SBATCH --array=1-48

######################
# Begin work section #
######################

$tissue=`sed "${SLURM_ARRAY_TASK_ID}q;d" tissuenames.txt`

echo "$tissue"

ml R
ml gcc

Rscript /work-zfs/rmccoy22/aseyedi2/neanderthal-sqtl/src/primary/R/NRich.R \
  "${tissue}" \
  "/work-zfs/rmccoy22/aseyedi2/GTExWGS_VCF/GTExWGS.AF.all.txt" \
  "/work-zfs/rmccoy22/aseyedi2/sqtl_permutation_backup/all_noms/varIDs/${tissue}_nom_varIDs.txt" \
  "/work-zfs/rmccoy22/aseyedi2/neanderthal-sqtl/analysis/SPRIME/sprime_calls.txt" \
  "/work-zfs/rmccoy22/aseyedi2/sqtl_permutation_backup/${tissue}_permutations.txt" \
  "/work-zfs/rmccoy22/aseyedi2/sqtl_permutation_backup/sqtl_nrich/${tissue}_enrichment.txt"
