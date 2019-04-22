#!/bin/bash
#SBATCH --partition=shared
#SBATCH --job-name=QTLTools-Loop
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

time { 
line=`sed "${SLURM_ARRAY_TASK_ID}q;d" GTExCovKey.csv`
full=$(echo $line | awk -F',' '{print $1}')
abb=$(echo $line | awk -F',' '{print $2}')
if grep "$abb" tissuesused.txt; then
   # Nominal Pass
   sbatch --wait --export=VCF=$VCF,pheno=$(echo $abb/$abb.pheno.bed.gz),tissue=$(echo $abb/$abb),covariates=$(echo $abb/$full.v7.covariates_output.txt) ${scripts}/sh/NomPass.sh

   for i in {1..100}; do
      cat $abb/${abb}_nominals_chunk_${i}.txt | gzip -c >> $abb/$abb.nominals.all.chunks.txt.gz
   done

   ls $abb/${abb}_* | sort -V >> $abb/${abb}_chunks_list.txt

   #Extract Neanderthal sequences
   cp ${data}/../analysis/SPRIME/sprime_calls.txt $abb
   sbatch --wait --export=listPath=$PWD/$abb,tissue=$(echo $abb),scripts=$scripts ${scripts}sh/NomPassExtractCall.sh
   rm $abb/sprime_calls.txt

   for i in {1..100}; do
      cat ${abb}/${abb}_nominals_chunk_${i}_out.txt | gzip -c >> ${abb}/$abb.nominals.all.chunks.NE_only.txt.gz
   done

   mkdir $abb/nominals; mv $abb/*_nominals_* $abb/nominals/

   #Call permuatation pass
   sbatch --wait --export=VCF=$VCF,pheno=$(echo $abb/$abb.pheno.bed.gz),tissue=$(echo $abb/$abb),covariates=$(echo $abb/$full.v7.covariates_output.txt) ${scripts}/sh/PermPass.sh

   for i in {1..100}; do
      cat $abb/${abb}_permutations_chunk_$i.txt | gzip -c >> $abb/${abb}.permutations_full.txt.gz
   done

   mkdir ${abb}/permutations; mv ${abb}/*_permutations_* ${abb}/permutations/

   Rscript ~/work/progs/QTLtools/script/runFDR_cis.R $abb/$abb.permutations_full.txt.gz 0.05 $abb/$abb.permutations_full_FDR

   sbatch --wait --export=VCF=$VCF,pheno=$(echo $abb/$abb.pheno.bed.gz),tissue=$(echo $abb/$abb),covariates=$(echo $abb/$full.v7.covariates_output.txt),permutations=$(echo $abb/$abb.permutations_full_FDR.thresholds.txt) ${scripts}/sh/CondPass.sh

   for i in {1..100}; do
      cat $abb/${abb}_conditionals_chunk_$i.txt | gzip -c >> $abb/${abb}.conditional_full.txt.gz
   done

   mkdir ${abb}/conditionals; mv ${abb}/*_conditionals_* ${abb}conditionals/

   mkdir /home-1/aseyedi2@jhu.edu/work/aseyedi2/sQTL/$abb

   Rscript ${scripts}/R/QQPlot-Viz.R /home-1/aseyedi2@jhu.edu/work/aseyedi2/sQTL/$abb $abb/$abb.nominals.all.chunks.NE_only.txt.gz $abb/$abb.permutations_full.txt.gz ${data}/../analysis/SPRIME/sprime_calls.txt
fi
echo "$full is done with QTL mapping"
}