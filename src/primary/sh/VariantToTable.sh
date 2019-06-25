#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=VariantToTable
#SBATCH --nodes=1
#SBATCH --time=0:30:0
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
#SBATCH --array=1-22

######################
# Begin work section #
######################

ml java
java -jar ~/work/progs/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar \
  -R Homo_sapiens_assembly19.fasta \
  -T VariantsToTable \
  -V GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz \
  -L ${SLURM_ARRAY_TASK_ID} \
  -F CHROM \
  -F POS \
  -F REF \
  -F ALT \
  -F AF \
  -o GTExWGS.AF.chr${SLURM_ARRAY_TASK_ID}.txt