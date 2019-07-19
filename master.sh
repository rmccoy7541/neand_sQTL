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
# This is the master script that prepares and submits all jobs for LeafCutter - all you need to do is put the GTEx sra's in the sra/ path#
#																																		 #
# Fill out the variables below with the appropriate full paths for each of the corresponding directories.                                #
# Also make sure to adjust the interactive session below for how long you think the whole pipeline will take. Err on the side of longer. #
#																																		 #
# It is also assumed that the GTEx VCF file for the whole genome has already been downloaded (and processed) and resides in ncbi/files/	 #
# 	(see Documentation for details)																										 #
#																																		 #
# It's also worth mentioning that though I've made variables for the pipeline, many of the shell scripts use absolute paths when calling #
# programs such as samtools, QTLtools, sra-tools etc. due to restraints in the HPC. I would recommend going through the shell scripts    #
# adapting them to your specific situation.																								 #
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
# input directory with sra files here
sra=$(echo /home-1/aseyedi2@jhu.edu/work/Ne_sQTL/sra/lung_skinEx_thy)
# leafcutter directory here
leafCutter=$(echo /scratch/groups/rmccoy22/aseyedi2/leafcutter)
# sprime
sprime=$(echo ${data}/../analysis/SPRIME/sprime_calls.txt)

## Step 1a - Conversion & Validation
################################################
# convert .sra to .bam files
cd $sra
# store all .sra names into text file for job array
ls *.sra >> sralist.txt

sra2BamNum=$(wc -l sralist.txt | awk '{print $1}')
# sra2bam, most computationally intensive step
sbatch --wait -a 1-$sra2BamNum --export=sraListPath=$PWD,homeDir=$homeDir ${scripts}/sh/sra2bam.sh
## samtools error check, remove broken bams
samtools quickcheck *bam 2> samtools_err_bam.txt

cat samtools_err_bam.txt | awk -F' ' '{print $1}' > failedbams.txt

for i in $(cat failedbams.txt)
do
   echo "$i is broken, removing now..." 
   rm $i
done

## Step 1b - Filtration
################################################

# list of bams to be filtered 
ls *.bam >> bamlist.txt
# bring bed file to current directory
cp ${data}/12-07-2018/GRCh37.bed $PWD
# filter unplaced contigs
bamNum=$(wc -l bamlist.txt | awk '{print $1}')

echo "Filtering unplaced contigs..."
sbatch --wait -a 1-$bamNum ${scripts}/sh/filter_bam.sh

samtools quickcheck *filt 2> samtools_err_filt.txt

cat samtools_err_filt.txt | awk -F' ' '{print $1}' > failedfilt.txt

for i in $(cat failedfilt.txt)
do
   echo "$i is broken, removing now..." >> log
   rm $i
done

ls *.filt >> filtlist.txt


## Step 2 - Intron Clustering
################################################
mkdir juncfiles
## IMPORTANT: changed leafcutter's bam2junc.sh to directly call
## the bin for samtools
filtNum=$(wc -l filtlist.txt | awk '{print $1}')

echo "Converting bam files to junc..."
sbatch --wait -a 1-$filtNum ${scripts}/sh/bam2junccall.sh
mv *.junc juncfiles/
cd juncfiles/
# strip junc files - STILL WITH RUN ID 'SRR######'
find -type f -name '*.sra.bam.filt.junc' | while read f; do mv "$f" "${f%.sra.bam.filt.junc}"; done

# put all of the renamed junc files in a text
ls SRR* > juncfiles.txt
# intron clustering
mkdir intronclustering/
echo "Intron clustering..."
# intron clustering
python $leafCutter/clustering/leafcutter_cluster.py -j juncfiles.txt -r intronclustering/ -m 50 -o Ne-sQTL -l 500000 # no option for parallelization
cd intronclustering/

## Step 3 - PCA calculation
################################################
# IMPORTANT: altered to remove 0 samples with NA's
echo "Preparing phenotype table..."
python $leafCutter/scripts/prepare_phenotype_table.py Ne-sQTL_perind.counts.gz -p 10 

# indexing and bedding
sh Ne-sQTL_perind.counts.gz_prepare.sh
# about the above line: you need to remove all of the index files and generate new ones once you convert the beds to a QTLtools compatible format

## Step 4 - VCF Preparation (optional, see doc for details)
################################################

## Step 5 - QTLtools Preparation
################################################
# prepare files for QTLtools
ls *qqnorm*.gz >> leafcutterphenotypes.txt 
# important: render these files compatible with QTLtools
echo "Making phenotype files QTLtools compatible..."
sbatch --wait ${scripts}/sh/QTLtools-Filter.sh
ls *.qtltools >> qtltools-input.txt
# generate the corresponding tbi files
rm Ne*tbi
for i in {1..22}; do tabix -p bed Ne-sQTL_perind.counts.gz.qqnorm_chr${i}.gz.qtltools; echo "Bedding chromosome $i"; done

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

# cats the nominal pass results
sbatch $scripts/primary/sh/CatNomPass.sh

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

# for i in *.txt; do
#    Rscript ${scripts}/R/count_sqtl.R #Something like that
# done

for i in $(ls *permutations.txt); do
   q=$(echo $i | cut -d _ -f 1)
   Rscript /scratch/groups/rmccoy22/progs/QTLtools/script/runFDR_cis.R $i 0.10 ${q}.results
done

for i in $(ls *results.significant.txt); do
  echo "Fixing $i to $i.bed..."
  cat $i | awk '{ print $9, $10-1, $11, $8, $1, $5 }' | tr " " "\t" | sort -k1,1 -k2,2n > $i.bed
done

for i in $(ls *_permutations.txt | sort -V); do echo $i | cut -d'_' -f 1; done > tissuenames.txt

for i in $(cat tissuenames.txt); do 
  echo "Submitting $i Enrichment Test"
  sbatch --export=tissue=$(echo $i) ${scripts}/sh/CallNRich.sh 
done
