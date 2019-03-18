#!/bin/bash

##########################################################################################################################################
# This is the master script that prepares and submits all jobs for LeafCutter								 #
# At the highest directory, it is assumed that the files are arranged as such:								 #
# Ne_sQTL/		aseyedi2/leafcutter/		aseyedi2/neand_sQTL/		master.sh					 #
#																	 #
# Don't ask why the git project name and the ncbi folder name (Ne_sQTL) are so similar. I messed up and am too afraid to fix it.	 #
#																	 #
# It is also assumed that the GTEx VCF file for the whole genome has already been downloaded and resides in ncbi/files/ 		 #
# 	(see Documentation for details)													 #
# Finally, please make sure that you have also downloaded the SRA files in ncbi/sra/ AND their corresponding SRA metadata 		 #
#														(SraRunTable.txt)	 #
# 	(again, please see Documentation for details)											 #
##########################################################################################################################################

# load modules
ml samtools
ml sra-tools
ml python/2.7-anaconda
ml bcftools
ml htslib
ml R
# ml qtltools - once MARCC guys install it
ml

# were are at ~ right now, which is /home-1/aseyedi2@jhu.edu
cd ~/work/

# the directory of master.sh
homeDir=$(echo ~/work/)
scripts=$(echo /home-1/aseyedi2@jhu.edu/work/aseyedi2/neand_sQTL/src/primary/)
data=$(echo /home-1/aseyedi2@jhu.edu/work/aseyedi2/neand_sQTL/data/)
ncbiFiles=$(echo /scratch/groups/rmccoy22/Ne_sQTL/files/)

# IF YOU ALREADY HAVE NON-BIALLELIC INDEXED VCF, uncomment below
VCF=$(echo /scratch/groups/rmccoy22/Ne_sQTL/files/GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz)

# input directory with sra files here
# Do this step for each separate batch of tissue samples you are trying to process concurrently
sra=$(echo /scratch/groups/rmccoy22/Ne_sQTL/sra/frontallobe_liver_muscle)

## Step 1a - Conversion & Validation
################################################
# convert .sra to .bam files
cd $sra
# store all .sra names into text file for job array
ls *.sra >> sralist.txt
# submit batch job, return stdout in $RES
sbatch --wait --export=sraListPath=$PWD,homeDir=$homeDir ${scripts}/sh/sra2bam.sh

samtools quickcheck *bam

### break here - deal accordingly with the bams

## Step 1b - Filtration
################################################

# list of bams to be filtered 
ls *.bam >> bamlist.txt
# bring bed file to current directory
cp ${data}/12-07-2018/GRCh37.bed $PWD
# filter unplaced contigs
sbatch --wait ${scripts}/sh/filter_bam.sh

## samtools quickcheck

ls *.filt >> filtlist.txt

## maybe inclue an if-statement after each sbatch that would catch any non-zero exit codes and abort the program

## Step 2 - Intron Clustering
################################################
mkdir juncfiles
## IMPORTANT: changed leafcutter's bam2junc.sh to directly call
## the bin for samtools
sbatch --wait ${scripts}/sh/bam2junccall.sh
mv *.junc juncfiles/
cd juncfiles/
# strip junc files - STILL WITH RUN ID 'SRR######'
find -type f -name '*.sra.bam.filt.junc' | while read f; do mv "$f" "${f%.sra.bam.filt.junc}"; done
# put all of the renamed junc files in a text

## make it v7 Analysis freeze juncs ONLY && GENOTYPED

ls SRR* >> juncfiles.txt
# intron clustering
mkdir intronclustering/
### call an interactive session with a good deal of memory for this step
python $homeDir/aseyedi2/leafcutter/clustering/leafcutter_cluster.py -j juncfiles.txt -r intronclustering/ -m 50 -o Ne-sQTL -l 500000 # no option for parallelization
cd intronclustering/

## Step 3 - PCA calculation
################################################
python $homeDir/aseyedi2/leafcutter/scripts/prepare_phenotype_table.py Ne-sQTL_perind.counts.gz -p 10 # works fine; also no option for parallelization
# indexing and bedding
ml htslib; sh Ne-sQTL_perind.counts.gz_prepare.sh
# about the above line: you need to remove all of the index files and generate new ones once you convert the beds to a QTLtools compatible format


############ Separate sh

## Step 4 - Genotype & Covariates Preparation
################################################
# filter the non-biallelic sites from genotype file using bcftools; see script for details
sbatch --wait ${scripts}/sh/bcf_tools.sh $homeDir
# index our friend with tabix
echo "Indexing our friend..."
ml htslib; tabix -p vcf GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz

############ Separate sh


## Step 5 - QTLtools Preparation
################################################


# prepare files for QTLtools
ls *qqnorm*.gz >> leafcutterphenotypes.txt 
# important: render these files compatible with QTLtools
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
for i in 1_*.txt; do echo $i | cut -d'_' -f 2| cut -d'.' -f 1 > tissuesused.txt; done

# moves each outputted file into its respective tissue folder
for i in *_*.txt; do echo $i | awk -F'[_.]' '{print $2}' | xargs -I '{}' mv $i '{}' ; done

# figure out how to make this part generalizable, but for now make a directory for each type of tissue processed

## figure out how to make this part generalizable too
ml htslib
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

ml bedtools

for line in $(cat tissuesused.txt)
do
   echo "Sorting $line.phen_fastqtl.bed to $line/$line.pheno.bed..."
   bedtools sort -header -i $line/$line.phen_fastqtl.bed > $line/$line.pheno.bed 
   echo "bgzipping $line/$line.pheno.bed..."
   bgzip $line/$line.pheno.bed
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

cp $data/Metadata/GTExCovKey.csv .

# Moves covariates to corresponding directory
for line in $(cat GTExCovKey.csv)
do
   full=$(echo $line | awk -F',' '{print $1}')
   abb=$(echo $line | awk -F',' '{print $2}')
   if grep "$abb" tissuesused.txt; then
      cp GTEx_Analysis_v7_eQTL_covariates/$full.v7.covariates.txt $abb
      Rscript ${scripts}/R/mergePCs.R Ne-sQTL_perind.counts.gz.PCs $abb/$full.v7.covariates.txt tissuetable/tissue_table.txt
   fi
done

for line in $(cat tissuesused.txt)
do

echo "Concatenating covariates..."

# this code is an absolute mess. I need to clean it up.

## Step 4 - Mapping sQTLs using QTLtools
################################################


##### Make this part useable for any tissue and not just whole blood
VCF=$(echo /scratch/groups/rmccoy22/Ne_sQTL/files/GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz)

#for loop for QTLtools nominals - Make this into a batch script -MAKE GENERALIZABLE

sbatch --wait --export=VCF=$ncbiFiles/GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz,pheno=$pheno ${scripts}/sh/NomPass.sh

# MAKE GENERALIZEABLE
cat WHLBLD_nominals_chunk_*.txt | gzip -c > nominals.all.chunks.txt.gz

ls WHLBLD_* | sort -V >> WHLBLD_chunks_list.txt

#Extract Neanderthal sequences
cp ${data}02-11-2019/tag_snps.neand.EUR.bed $PWD
sbatch --export=listPath=$PWD,tissue=$(echo BRNCHA),scripts=$scripts ${scripts}sh/NomPassExtractCall.sh

cat BRNCHA_nominals_chunk_*_out.txt | gzip -c > BRNCHA.nominals.all.chunks.NE_only.txt.gz

mkdir nominals; mv *_nominals_* nominals/

#Call permuatation pass
sbatch --wait --export=VCF=$VCF,pheno=$pheno,tissue=$(echo BRNCHA),covariate=$(echo Brain_Cerebellum.v7.covariates_output.txt) ${scripts}/sh/PermPass.sh

cat TESTIS_permutations_chunk_*.txt | gzip -c > TESTIS.permutations_full.txt.gz

mkdir permutations; mv *_permutations_* permutations/

ml gcc

Rscript ~/work/progs/QTLtools/script/runFDR_cis.R BRNCHA.permutations_full.txt.gz 0.05 BRNCHA.permutations_full_FDR

sbatch --wait --export=VCF=$VCF,pheno=$pheno,tissue=$(echo TESTIS),covariates=$(echo Testis.v7.covariates_output.txt),permutations=$(echo TESTIS.permutations_full_FDR.thresholds.txt) ${scripts}/sh/CondPass.sh

sbatch --wait --export=VCF=$VCF,pheno=$pheno,tissue=$(echo TESTIS),covariates=$(echo Brain_Cerebellum.v7.covariates_output.txt),permutations=$(echo TESTIS.permutations_full_FDR.thresholds.txt) ${scripts}/sh/CondPass.sh

cat WHLBLD_conditionals_chunk_*.txt | gzip -c > conditional_full.txt.gz

mkdir conditionals; mv *_conditionals_* conditionals/

Rscript ${scripts}/R/QQPlot-Viz.R /home-1/aseyedi2@jhu.edu/work/aseyedi2/sQTL/WHLBLD WHLBLD.nominals.all.chunks.NE_only.txt.gz WHLBLD.permutations_full.txt.gz

echo "done"
