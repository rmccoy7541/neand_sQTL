#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=NomPass
#SBATCH --nodes=1
#SBATCH --time=4:0:0
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --array=1-100

######################
# Begin work section #
######################

# this project's scripts dir
scripts=$(echo /home-1/aseyedi2@jhu.edu/work/aseyedi2/neand_sQTL/src/primary/)
# data dir
data=$(echo /home-1/aseyedi2@jhu.edu/work/aseyedi2/neand_sQTL/data/)
# ncbi/files/
ncbiFiles=$(echo /scratch/groups/rmccoy22/Ne_sQTL/files/)
# IF YOU ALREADY HAVE NON-BIALLELIC INDEXED VCF
VCF=$(echo /scratch/groups/rmccoy22/Ne_sQTL/files/GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz)
# input directory with junc files here
junc=$(echo '/scratch/groups/rmccoy22/Ne_sQTL/sra/sqtl_junc')
# leafcutter directory here
leafCutter=$(echo /scratch/groups/rmccoy22/aseyedi2/leafcutter)

/scratch/groups/rmccoy22/progs/QTLtools/QTLtools_1.1_Ubuntu14.04_x86_64 cis \
    --vcf  $VCF \
    --bed "$pheno" \
    --cov  "$covariates" \
    --nominal 1  \
    --chunk $SLURM_ARRAY_TASK_ID 100 \
    --out "${tissue}_nominals_chunk_${SLURM_ARRAY_TASK_ID}.txt"

scripts=$(echo /home-1/aseyedi2@jhu.edu/work/aseyedi2/neand_sQTL/src/primary/)

sbatch -a $SLURM_ARRAY_TASK_ID --export=listPath=$PWD/$abb,tissue=$(echo $abb),scripts=$scripts ${scripts}sh/NomPassExtractCall.sh

#Call permuatation pass
sbatch -a $SLURM_ARRAY_TASK_ID --export=VCF=$VCF,pheno=$(echo $abb/$abb.pheno.bed.gz),tissue=$(echo $abb/$abb),covariates=$(echo $abb/$full.v7.covariates_output.txt),abb=$abb,full=$full ${scripts}/sh/PermPass.sh