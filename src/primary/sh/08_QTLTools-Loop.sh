#!/bin/bash
#SBATCH --partition=shared
#SBATCH --job-name=QTLTools-Loop
#SBATCH --time=0:1:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

ml htslib
ml R
ml gcc


line=`sed "${SLURM_ARRAY_TASK_ID}q;d" GTExCovKey.csv`

full=$(echo $line | awk -F',' '{print $1}')
abb=$(echo $line | awk -F',' '{print $2}')

cp $sprime $abb
# This next script does nom pass, then calls perm pass AND nom pass extract
sbatch --export=VCF=$VCF,pheno=$(echo $abb/$abb.pheno.bed.gz),tissue=$(echo $abb/$abb),covariates=$(echo $abb/$full.v7.covariates_output.txt),abb=$abb,full=$full,worDir=$PWD,scripts=$scripts,data=$data,sprime=$sprime \
  ${scripts}/sh/09_NomPass.sh