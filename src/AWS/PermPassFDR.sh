#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=PermPassFDR_cis
#SBATCH --nodes=1
#SBATCH --time=2:0:0
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3

######################
# Begin work section #
######################

Rscript ~/work/progs/QTLtools/script/runFDR_cis.R $abb/$abb.permutations_full.txt.gz 0.05 $abb/$abb.permutations_full_FDR