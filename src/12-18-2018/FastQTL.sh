#!/bin/bash
#SBATCH --job-name=FastQTL
#SBATCH --time=2:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-22


line=`sed "${SLURM_ARRAY_TASK_ID}q;d" FastQTLPhenotypes.txt`
chr = $(expr "/$line" : '.*\([^/.]\{4\}\)\.[^/.]*$')
for i in {1..10}
do
	fastQTL -V GTExWholeGenomeSequenceGenotypeMatrixBiallelicOnly.vcf.gz -B $line -C testNE_sQTL_perind.counts.gz.PCs -O ${chr}_${i}_FastQTLresult.txt -L ${chr}.log --chunk ${i} 10
done
