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

ml R
ml gcc

Rscript ~/work/progs/QTLtools/script/runFDR_cis.R $abb/$abb.permutations_full.txt.gz 0.05 $abb/$abb.permutations_full_FDR

sbatch -a 1-100 --export=VCF=$VCF,pheno=$(echo $abb/$abb.pheno.bed.gz),tissue=$(echo $abb/$abb),covariates=$(echo $abb/$full.v7.covariates_output.txt),permutations=$(echo $abb/$abb.permutations_full_FDR.thresholds.txt),abb=$abb,full=$full ${scripts}/sh/CondPass.sh   
