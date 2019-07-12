#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=NRich.sh
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24

######################
# Begin work section #
######################

ml R

date

Rscript /work-zfs/rmccoy22/aseyedi2/neanderthal-sqtl/src/primary/R/NRich.R \
  "/work-zfs/rmccoy22/aseyedi2/GTExWGS_VCF/GTExWGS.AF.all.txt" \
  "/work-zfs/rmccoy22/rmccoy22/sqtl/intron_clustering/tryagain/" \
  "/work-zfs/rmccoy22/aseyedi2/neanderthal-sqtl/analysis/SPRIME/sprime_calls.txt" \
  "/work-zfs/rmccoy22/aseyedi2/GTExWGS_VCF/${tissue}_permutations.txt" \
  "/work-zfs/rmccoy22/aseyedi2/sqtl_permutation_backup/sqtl_nrich/${tissue}_enrichment.txt"

date