#!/bin/bash
#SBATCH --job-name=index_vcf
#SBATCH --time=6:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=6

# index using tabix
ml htslib; tabix -p vcf $outdir/GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz