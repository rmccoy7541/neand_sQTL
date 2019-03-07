#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=NomPass
#SBATCH --nodes=1
#SBATCH --time=8:0:0
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --array=1-100

######################
# Begin work section #
######################

/scratch/groups/rmccoy22/progs/QTLtools/QTLtools_1.1_Ubuntu14.04_x86_64 cis \
	--vcf  $VCF \
	--bed "$pheno" \
	--cov  "Brain_Cerebellum.v7.covariates_output.txt" \
	--nominal 1  \
	--chunk $i 100 \
	--out "BRNCHA_nominals_chunk_${i}.txt"
done
