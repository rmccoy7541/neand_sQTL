#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=VariantToTable
#SBATCH --nodes=1
#SBATCH --time=6:0:0
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --array=1-22

######################
# Begin work section #
######################

# extracts allele frequences from the VCF
ml java
java -jar GenomeAnalysisTK \
  -R Homo_sapiens_assembly19.fasta \
  -T VariantsToTable \
  -V $VCF \
  -L ${SLURM_ARRAY_TASK_ID} \
  -F CHROM \
  -F POS \
  -F REF \
  -F ALT \
  -F AF \
  --showFiltered \
  -o GTExWGS.AF.chr${SLURM_ARRAY_TASK_ID}.txt
