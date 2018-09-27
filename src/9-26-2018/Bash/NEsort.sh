#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=sortNE
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=end
#SBATCH --mail-user=aseyedi2@jhu.edu
#SBATCH --array=3

######################
# Begin work section #
######################

awk '{print $1"_"$3"\t"$4"\t"$5}' tag_snps.neand.EAS.bed | sort -k1,1 > "tag_snps.neand.EAS.sorted.txt"
awk '{print $1"_"$3"\t"$4"\t"$5}' tag_snps.neand.EUR.bed | sort -k1,1 > "tag_snps.neand.EUR.sorted.txt"
awk '{print $1"_"$3"\t"$4"\t"$5}' tag_snps.neand.SAS.bed | sort -k1,1 > "tag_snps.neand.SAS.sorted.txt"