#!/bin/bash
#SBATCH --partition=unlimited
#SBATCH --job-name=master-shell
#SBATCH --time=14-0:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
######################
# Begin work section #
##########################################################################################################################################
#  ███╗   ██╗███████╗ █████╗ ███╗   ██╗██████╗ ███████╗██████╗ ████████╗██╗  ██╗ █████╗ ██╗         ███████╗ ██████╗ ████████╗██╗     
#  ████╗  ██║██╔════╝██╔══██╗████╗  ██║██╔══██╗██╔════╝██╔══██╗╚══██╔══╝██║  ██║██╔══██╗██║         ██╔════╝██╔═══██╗╚══██╔══╝██║     
#  ██╔██╗ ██║█████╗  ███████║██╔██╗ ██║██║  ██║█████╗  ██████╔╝   ██║   ███████║███████║██║         ███████╗██║   ██║   ██║   ██║     
#  ██║╚██╗██║██╔══╝  ██╔══██║██║╚██╗██║██║  ██║██╔══╝  ██╔══██╗   ██║   ██╔══██║██╔══██║██║         ╚════██║██║▄▄ ██║   ██║   ██║     
#  ██║ ╚████║███████╗██║  ██║██║ ╚████║██████╔╝███████╗██║  ██║   ██║   ██║  ██║██║  ██║███████╗    ███████║╚██████╔╝   ██║   ███████╗
#  ╚═╝  ╚═══╝╚══════╝╚═╝  ╚═╝╚═╝  ╚═══╝╚═════╝ ╚══════╝╚═╝  ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝    ╚══════╝ ╚══▀▀═╝    ╚═╝   ╚══════╝
##########################################################################################################################################
# Fill out the variables below with the appropriate full paths for each of the corresponding directories.                                #
# Also make sure to adjust the interactive session below for how long you think the whole pipeline will take. Err on the side of longer. #
#																																		                                                                     #
# It is also assumed that the GTEx VCF file for the whole genome has already been downloaded (and processed) and resides in ncbi/files/	 #
# 	(see Documentation for details)																										                                                   #
#																																		                                                                     #
# It's also worth mentioning that though I've made variables for the pipeline, many of the shell scripts use absolute paths when calling #
# programs such as samtools, QTLtools, sra-tools etc. due to restraints in the HPC. I would recommend going through the shell scripts    #
# adapting them to your specific situation.																								                                               #
##########################################################################################################################################

# adjust for the interactive session time length you might need here

# load modules
ml bedtools
ml samtools
ml sra-tools
ml python/2.7-anaconda
ml bcftools
ml htslib
ml R
ml gcc
ml

#### Fill the Directories appropriately here
# top-level directory, above ncbi/
homeDir=$(echo ~/work/)
# this project's scripts dir
scripts=$(echo /work-zfs/rmccoy22/aseyedi2/neanderthal-sqtl/src/primary/)
# data dir
data=$(echo /work-zfs/rmccoy22/aseyedi2/neanderthal-sqtl/data/)
# ncbi/files/
ncbiFiles=$(echo /scratch/groups/rmccoy22/Ne_sQTL/files/)
# IF YOU ALREADY HAVE NON-BIALLELIC INDEXED VCF
VCF=$(echo /scratch/groups/rmccoy22/Ne_sQTL/files/GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz)
# leafcutter directory here
leafCutter=$(echo /scratch/groups/rmccoy22/aseyedi2/leafcutter)
# sprime
sprime=$(echo ${data}/../analysis/SPRIME/sprime_calls.txt)


mkdir intronclustering/

sbatch --wait ${scripts}/../AWS/junc_cluster.sh

cd intronclustering/

echo "Preparing phenotype table..."
sbatch --wait ${scripts}/../AWS/prepare_phen_table.sh $leafCutter

## Step 5 - QTLtools Preparation
################################################
# prepare files for QTLtools
ls *qqnorm* >> leafcutterphenotypes.txt 
# important: render these files compatible with QTLtools
echo "Making phenotype files QTLtools compatible..."
sbatch --wait ${scripts}/sh/QTLtools-Filter.sh
ls *.qtltools >> qtltools-input.txt

# generate the corresponding tbi files
interact -p shared -t 04:0:0 -c 3
for i in {1..22}; do tabix -p bed Ne-sQTL_perind.counts.gz.qqnorm_chr${i}.gz.qtltools; echo "Bedding chromosome $i"; done
exit

cp ${data}/01-22-2019/GTExTissueKey.csv $PWD
# get the tissue sites for each corresonding sra file
Rscript ${scripts}/R/sraTissueExtract.R ${data}/Metadata/SraRunTable.txt GTExTissueKey.csv

# submit each LF phenotype file to sraNameChangeSort as command line variable as well as tissue_table.txt
for phen in *qqnorm*.gz.qtltools; do Rscript ${scripts}/R/sraNameChangeSort.R $phen tissue_table.txt ; done
cat tissue_table.txt | cut -f3 | awk '{if(NR>1)print}' |  awk '!seen[$0]++' > tissuenames.txt

mkdir tissuetable/
mv tissue_table.txt tissuetable/
# make directories for each type of tissue
for i in 1_*.txt; do echo $i | cut -d'_' -f 2| cut -d'.' -f 1 | xargs mkdir; done

# save tissue types in a file
for i in 1_*.txt; do echo $i | cut -d'_' -f 2| cut -d'.' -f 1 >> tissuesused.txt; done

# moves each outputted file into its respective tissue folder
for i in *_*.txt; do echo $i | awk -F'[_.]' '{print $2}' | xargs -I '{}' mv $i '{}' ; done

## Concatting the phenotype files
for line in $(cat tissuesused.txt)
do
   head -1 $line/1_$line.txt > $line/$line.phen_fastqtl.bed
   echo "Concatenating $line phenotypes..."
   for file in $(ls $line/*_*.txt)
   do
      echo "Adding $file..."
      cat $file | sed -e1,1d >> $line/$line.phen_fastqtl.bed
   done
done

for line in $(cat tissuesused.txt)
do
   echo "Sorting $line.phen_fastqtl.bed to $line/$line.pheno.bed..."
   bedtools sort -header -i $line/$line.phen_fastqtl.bed > $line/$line.pheno.bed 
   echo "bgzipping $line/$line.pheno.bed..."
   bgzip -f $line/$line.pheno.bed
   #figure out where tabix outputs
   echo "Indexing $line/$line.pheno.bed.gz..."
   tabix -p bed $line/$line.pheno.bed.gz
done

for line in $(cat tissuesused.txt)
do
   mkdir $line/sepfiles
   mv $line/*_${line}.txt $line/sepfiles/
done

# download genotype covariates
wget https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL_covariates.tar.gz

tar -xvf GTEx_Analysis_v7_eQTL_covariates.tar.gz

cp $data/Metadata/GTExCovKey.csv $PWD

# Moves covariates to corresponding directory
for line in $(cat GTExCovKey.csv)
do
   full=$(echo $line | awk -F',' '{print $1}')
   abb=$(echo $line | awk -F',' '{print $2}')
   if grep "$abb" tissuesused.txt; then
      cp GTEx_Analysis_v7_eQTL_covariates/$full.v7.covariates.txt $abb
      Rscript ${scripts}/R/mergePCs.R Ne-sQTL_perind.counts.gz.PCs $abb/$full.v7.covariates.txt tissuetable/tissue_table.txt
      mv $full.v7.covariates_output.txt $abb
   fi
done

## Step 4 - Mapping sQTLs using QTLtools
################################################
numTissues=$(wc -l GTExCovKey.csv)

# Will take at least 3 weeks lol
sbatch --wait -a 2-$numTissues ${scripts}/sh/QTLTools-Loop.sh

mkdir /scratch/groups/rmccoy22/rmccoy22/sqtl_permutation_backup
mkdir all_noms
mkdir sqtl_nrich
cd /scratch/groups/rmccoy22/rmccoy22/sqtl_permutation_backup
cp /work-zfs/rmccoy22/rmccoy22/sqtl/intron_clustering/tryagain/*/*permutation* .
cp /work-zfs/rmccoy22/rmccoy22/sqtl/intron_clustering/tryagain/*/*pheno* .
cp /work-zfs/rmccoy22/rmccoy22/sqtl/intron_clustering/tryagain/tissuenames.txt .
cp ${data}/../analysis/SPRIME/sprime_calls.txt .

# concatenates permutation pass results
for TISSUE in ADPSBQ ADPVSC ADRNLG ARTACRN ARTAORT ARTTBL BREAST BRNACC BRNAMY BRNCDT BRNCHA BRNCHB BRNCTXA BRNCTXB BRNHPP BRNHPT BRNNCC BRNPTM BRNSNG BRNSPC CLNSGM CLNTRN ESPGEJ ESPMCS ESPMSL FIBRLBLS HRTAA HRTLV LCL LIVER LUNG MSCLSK NERVET OVARY PNCREAS PRSTTE PTTARY SKINNS SKINS SLVRYG SNTTRM SPLEEN STMACH TESTIS THYROID UTERUS VAGINA WHLBLD
do
  mkdir ${TISSUE}
  mv ${TISSUE}*permutations* ${TISSUE}
  for i in {1..100}
  do
    echo "Catting ${TISSUE}_permutations_chunk_${i}.txt..."
    cat ${TISSUE}/${TISSUE}_permutations_chunk_${i}.txt >> ${TISSUE}_permutations.txt
  done
done

# getting tissue names
for i in $(ls *_permutations.txt | sort -V); do echo $i | cut -d'_' -f 1; done > tissuenames.txt

sbatch $scripts/sh/VariantToTable.sh

# concat the GTEx AF VCF chunks
cd /work-zfs/rmccoy22/aseyedi2/GTExWGS_VCF
cat GTExWGS.AF.chr1.txt > GTExWGS.AF.all.txt
for i in {2..22}; do tail -n +2 GTExWGS.AF.chr${i}.txt >> GTExWGS.AF.all.txt ; done
cd -

cd ~/data/aseyedi2/sqtl_permutation_backup/all_noms/varIDs/
sbatch --export=seed=$(echo "123"),M=$(echo "100") ${scripts}/sh/CallNRich.sh
