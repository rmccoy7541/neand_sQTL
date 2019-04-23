#!/bin/bash
#SBATCH --partition=shared
#SBATCH --job-name=QTLTools-Master
#SBATCH --time=12:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4

ml htslib
ml R
ml gcc
ml

# top-level directory, above ncbi/
homeDir=$(echo ~/work/)
# this project's scripts dir
scripts=$(echo /home-1/aseyedi2@jhu.edu/work/aseyedi2/neand_sQTL/src/primary/)
# data dir
data=$(echo /home-1/aseyedi2@jhu.edu/work/aseyedi2/neand_sQTL/data/)
# ncbi/files/
ncbiFiles=$(echo /scratch/groups/rmccoy22/Ne_sQTL/files/)
# IF YOU ALREADY HAVE NON-BIALLELIC INDEXED VCF
VCF=$(echo /scratch/groups/rmccoy22/Ne_sQTL/files/GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz)
# input directory with sra files here
sra=$(echo /home-1/aseyedi2@jhu.edu/work/Ne_sQTL/sra/lung_skinEx_thy)
# leafcutter directory here
leafCutter=$(echo /scratch/groups/rmccoy22/aseyedi2/leafcutter)


line=`sed "${SLURM_ARRAY_TASK_ID}q;d" GTExCovKey.csv`
full=$(echo $line | awk -F',' '{print $1}')
abb=$(echo $line | awk -F',' '{print $2}')
if grep "$abb" tissuesused.txt; then

   cp ${data}/../analysis/SPRIME/sprime_calls.txt $abb
   # This next script does nom pass, then calls perm pass AND nom pass extract, and then perm pass calls cond pass once finished, which cats everything or whatever.
   sbatch --wait --export=VCF=$VCF,pheno=$(echo $abb/$abb.pheno.bed.gz),tissue=$(echo $abb/$abb),covariates=$(echo $abb/$full.v7.covariates_output.txt) ${scripts}/sh/NomPass.sh

   


   mkdir ${abb}/conditionals; mv ${abb}/*_conditionals_* ${abb}conditionals/

   mkdir /home-1/aseyedi2@jhu.edu/work/aseyedi2/sQTL/$abb

   sbatch --wait --export=tissue=$(echo $abb) ${scripts}/../AWS/QTLtools-int.sh

   Rscript ${scripts}/R/QQPlot-Viz.R /home-1/aseyedi2@jhu.edu/work/aseyedi2/sQTL/$abb $abb/$abb.nominals.all.chunks.NE_only.txt.gz $abb/$abb.permutations_full.txt.gz ${data}/../analysis/SPRIME/sprime_calls.txt
   echo "$full is done with QTL mapping"
else
   echo "$full was not included in the analysis"
fi