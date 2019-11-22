#!/bin/bash -l
#SBATCH --job-name=tablemerge
#SBATCH -N 1
#SBATCH --time=15:00:00
#SBATCH --partition=lrgmem
#SBATCH --mem=400G
#SBATCH --ntasks-per-node=48
#SBATCH --mail-type=end
#SBATCH --mail-user=syan11@jhu.edu
#SBATCH --array=1-49

cd /scratch/groups/rmccoy22/syan11/sQTL_results/gtex_v8/merge

ml R/3.5.1; ml gcc

TISSUE=`sed "${SLURM_ARRAY_TASK_ID}q;d" tissues.txt`
INTRONS_FILE=../splitIC/${TISSUE}_intronCounts.txt
#INTRONS_FILE=../intron_counts/GTEx_v8_junctions_nohead.gct.gz
SQTL_FILE=/scratch/groups/rmccoy22/aseyedi2/sQTLv8/data/GTEx_Analysis_v8_sQTL/${TISSUE}.v8.sqtl_signifpairs.txt.gz
Rscript table_merge_SY.R $INTRONS_FILE $SQTL_FILE $TISSUE
