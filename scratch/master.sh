#!/bin/bash

##########################################################################################################################################
# This is the master script that prepares and submits all jobs for LeafCutter								 #
# At the highest directory, it is assumed that the files are arranged as such:								 #
# ncbi/		leafcutter/		Ne-sQTL/		master.sh								 #
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
# the directory of master.sh
homeDir = $(pwd -P)

## Step 1 - Conversion
################################################
# convert .sra to .bam files
cd ncbi/sra
# store all .sra names into text file for job array
ls *.sra >> sralist.txt
# submit batch job, return stdout in $RES
#jid1=$(sbatch --wait ${homeDir}/ncbi/src/12-14-2018/sra2bam.sh)
jid1=$(sbatch --wait ${homeDir}/Ne-sQTL/src/12-14-2018/sra2bam.sh)
# list of bams to be filtered 
ls *.bam >> bamlist.txt
# bring bed file to current directory
cp ${homeDir}/Ne-sQTL/data/12-07-2018/GRCh37.bed $PWD
# filter unplaced contigs
jid2=$(sbatch --wait --dependency=afterok:${jid1##* } ${homeDir}/Ne-sQTL/src/12-14-2018/filter_bam.sh)
ls *.filt >> filtlist.txt
# No longer renaming SRAs until after leafcutter
## maybe inclue an if-statement after each sbatch that would catch any non-zero exit codes and abort the program

## Step 2 - Intron Clustering
################################################
mkdir juncfiles
sbatch --wait --dependency=afterok:${jid2##* } ${homeDir}/Ne-sQTL/src/12-10-2018/bam2junccall.sh
mv *.junc juncfiles/
cd juncfiles/
# strip junc files - STILL WITH RUN ID 'SRR######'
find -type f -name '*.sra.bam.filt.junc' | while read f; do mv "$f" "${f%.sra.bam.filt.junc}"; done
# put all of the renamed junc files in a text
ls SRR* >> juncfiles.txt
# intron clustering
mkdir intronclustering/
### call an interactive session with a good deal of memory for this step
python $homeDir/leafcutter/clustering/leafcutter_cluster.py -j juncfiles.txt -r intronclustering/ -m 50 -o Ne-sQTL -l 500000 # no option for parallelization
cd intronclustering/

## Step 3 - PCA calculation
################################################
python $homeDir/leafcutter/scripts/prepare_phenotype_table.py Ne-sQTL_perind.counts.gz -p 10 # works fine; also no option for parallelization
# indexing and bedding
ml htslib; sh Ne-sQTL_perind.counts.gz_prepare.sh
# about the above line: you need to remove all of the index files and generate new ones once you convert the beds to a QTLtools compatible format


## Step 4 - Genotype & Covariates Preparation
################################################
# filter the non-biallelic sites from genotype file using bcftools; see script for details
sbatch --wait ${homeDir}/Ne-sQTL/src/12-18-2018/bcf_tools.sh $homeDir
# index our friend with tabix
echo "Indexing our friend..."
tabix -p vcf GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz
# prepare files for QTLtools
ls *qqnorm*.gz >> leafcutterphenotypes.txt 
# important: render these files compatible with QTLtools
sbatch --wait ${homeDir}/Ne-sQTL/src/12-18-2018/QTLtools-Filter.sh
ls *.qtltools >> qtltools-input.txt
# generate the corresponding tbi files
rm Ne*tbi
for i in {1..22}; do tabix -p bed Ne-sQTL_perind.counts.gz.qqnorm_chr${i}.gz.qtltools; echo "Bedding chromosome $i"; done
# download genotype covariates
wget https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL_covariates.tar.gz
tar -xzf GTEx_Analysis_v7_eQTL_covariates.tar.gz

## Step 5 - QTLtools Preparation
################################################

mv ${homeDir}/Ne-sQTL/data/01-22-2019/GTExTissueKey.txt $PWD
# get the tissue sites for each corresonding sra file
Rscript --vanilla ${homeDir}/Ne-sQTL/src/01-09-2019/sraTissueExtract.R ${homeDir}/Ne-sQTL/data/SraRunTable.txt GTExTissueKey.csv


# submit each LF phenotype file to sraNameChangeSort as command line variable as well as tissue_table.txt
for phen in *qqnorm*.gz.qtltools; do Rscript --vanilla ${homeDir}/Ne-sQTL/src/01-15-2019/sraNameChangeSort.R $phen tissue_table.txt ; done
rm *Leukemia*

#this code is problematic - I can't have any whitespaces in filenames, yet the that's all they have in tissue_table
# getting the tissue names from metadata - 48 tissues in all
cat tissue_table.txt | cut -f3 | awk '{if(NR>1)print}' |  awk '!seen[$0]++' > tissuenames.txt
mkdir cattissues

# concatenate all tissues from chr 1-22 such that all chromosomes are found in one file per tissue
for i in {1..48}; do line=`sed "${i}q;d" tissuenames.txt`; echo "Concatenating $line..."; for q in {1..22}; do echo "Chr $q..."; awk 'FNR==1 && NR!=1{next;}{print}' "${q}_${line}.txt" > "${line}".txt; rm "${q}_${line}.txt" done; done #formerly cat "${q}_${line}.txt" >> ... let's see if this works ; did this to prevent writing the header multiple times

mv 1-22_* cattissues/

# Parts of this next line may appear redundant. That's okay. it just creates a directory for each concatenated tissue phenotype file and moves it in there.
for file in cattissues/*; do i=$(echo $file); q=$(basename "$i"); filename="${q%.*}"; filename=$(echo ${filename#*_}); mkdir tissues/"${filename}"; mv "$file" tissues/"${filename}"; done
# at this point, you want to pass each tissue PC file and the leafcutter PC file as command-line arguments into an R script that concatenates the PCs by GTEX ID
for tissue in GTEx_Analysis_v7_eQTL_covariates/*; do Rscript --vanilla $homeDir/src/12-21-2018/mergePCs.R Ne-sQTL_perind.counts.gz.PCs ${tissue}; echo "Concatenating ${tissue}"; done
mv *covariates_* tissues/
cd tissues/

# Moving covariates to corresponding directories
################################################
#for i in *.txt;
#do
#  string=$i
#  until_dot=${string%%.*}
#  echo "$until_dot"
#
#  # Replace all non-alphabetic characters by the glob *
#  glob_pattern=${until_dot//[^[:alpha:]]/*}
#  echo "$glob_pattern"
#
#  # Use nullglob to have non matching glob expand to nothing
#  shopt -s nullglob
#  # DO NOT USE QUOTES IN THE FOLLOWING EXPANSION:
#  # the variable is actually a glob!
#  # Could also do dirs=( $glob_pattern*/ ) to check if directory
#  dirs=( $glob_pattern/ )
#
#  # Now check how many matches there are:
#  if ((${#dirs[@]} == 0)); then
#      echo >&2 "No matches for $glob_pattern"
#  elif ((${#dirs[@]} > 1)); then
#      echo >&2 "More than one matches for $glob_pattern: ${dirs[@]}"
#  else
#      echo "All good"
#      # Remove the echo to actually perform the move
#      echo mv "$string" "${dirs[0]}"
#  fi;
#done


## This doesn't work as it should - just going to do it manually


## Step 4 - Mapping sQTLs using QTLtools
################################################

find . -name "[1-22]*gz" -print0 | xargs -0 ls > phenpaths.txt


phenpaths=$(sed 's/ /\\ /g' <<< cat phenpaths.txt)

IFS=$'\n'       # make newlines the only separator

for tissue in $(cat ./phenpaths.txt)    
do
    phenpaths=$(sed 's/ /\\ /g' $tissue); echo $phenpaths 
done

# gets rid of duplicate headers, removes row numbers
for tissue in $(cat ./phenpaths.txt)    
do
    awk '{ if($0 != header) { print; } if(header == "") { header=$0; } }' "$tissue" | sed 's/^[0-9][0-9]*[ \t]*//' | awk '$3!=""' > "$tissue"
done

for tissue in $(cat ./phenpaths.txt)    
do
    bgzip -f "${tissue}"
done

#same as before, but with .gz
find . -name "1-22*" -print0 | xargs -0 ls > phenpaths.txt

for tissue in $(cat ./phenpaths.txt)    
do
    tabix -p bed $tissue
done

for tissue in $(cat ./phenpaths.txt)    
do
    for chunk in {1..20}; do sbatch --export=tissue=$tissue,chunk=$chunk QTLtools-NomPass.sh; done
done



