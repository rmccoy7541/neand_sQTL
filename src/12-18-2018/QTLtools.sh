#!/bin/bash
#SBATCH --job-name=QTLtools
#SBATCH --time=2:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-22

ml qtltools
# get file;
line=`sed "${SLURM_ARRAY_TASK_ID}q;d" QTLtools-input.txt`
# chromosome number
chr="chr${line//[!0-9]/}"
# FastQTL command
QTLtools -V GTExWholeGenomeSequenceGenotypeMatrixBiallelicOnly.vcf.gz -B $line -C testNE_sQTL_perind.counts.gz.PCs -O "${chr}_QTLtools-result.txt" -L "$chr.log" --chunk 1 1

