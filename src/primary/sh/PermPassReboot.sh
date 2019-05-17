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

# this project's scripts dir
scripts=$(echo /work-zfs/rmccoy22/aseyedi2/neanderthal-sqtl/src/primary)
# data dir
data=$(echo /work-zfs/rmccoy22/aseyedi2/neanderthal-sqtl/data)
# IF YOU ALREADY HAVE NON-BIALLELIC INDEXED VCF
VCF=$(echo /scratch/groups/rmccoy22/Ne_sQTL/files/GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz)
# sprime_calls
sprime=$(echo /work-zfs/rmccoy22/aseyedi2/neanderthal-sqtl/analysis/SPRIME/sprime_calls.txt)
# /work-zfs/rmccoy22/rmccoy22/sqtl/intron_clustering/tryagain

line=`sed "${SLURM_ARRAY_TASK_ID}q;d" GTExCovKey.csv`

full=$(echo $line | awk -F',' '{print $1}')
abb=$(echo $line | awk -F',' '{print $2}')


cp $sprime $abb
# This next script does nom pass, then calls perm pass AND nom pass extract, and then perm pass calls cond pass once finished, which cats everything or whatever.
#sbatch --export=VCF=$VCF,pheno=$(echo $abb/$abb.pheno.bed.gz),tissue=$(echo $abb/$abb),covariates=$(echo $abb/$full.v7.covariates_output.txt),abb=$abb,full=$full,worDir=$PWD ${scripts}/sh/NomPass.sh


# figure out how to implement these next two scripts.
sbatch -a 1-100 --export=VCF=$VCF,pheno=$(echo $abb/$abb.pheno.bed.gz),tissue=$(echo $abb/$abb),covariates=$(echo $abb/$full.v7.covariates_output.txt),permutations=$(echo $abb/$abb.permutations_full_FDR.thresholds.txt),abb=$abb,full=$full ${scripts}/sh/PermPass.sh

# for line in $(cat GTExCovKey.csv); do
	# full=$(echo $line | awk -F',' '{print $1}')
	# abb=$(echo $line | awk -F',' '{print $2}')
	# if grep $abb "failcat"; then
		# sbatch --export=VCF=$VCF,pheno=$(echo $abb/$abb.pheno.bed.gz),tissue=$(echo $abb/$abb),covariates=$(echo $abb/$full.v7.covariates_output.txt),permutations=$(echo $abb/$abb.permutations_full_FDR.thresholds.txt),abb=$abb,full=$full $scripts/../AWS/QTLtools-int.sh
	# fi
# done	