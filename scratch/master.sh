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
mv aseyedi2/neand_sQTL/master.sh $PWD

# the directory of master.sh
homeDir=$(pwd -P)
scripts=$("/home-1/aseyedi2@jhu.edu/work/aseyedi2/neand_sQTL/src/Primary/")

## Step 1 - Conversion
################################################
# convert .sra to .bam files
cd Ne_sQTL/sra
# store all .sra names into text file for job array
ls *.sra >> sralist.txt
# submit batch job, return stdout in $RES
sbatch --wait --export=sraListPath=$PWD ${scripts}/sh/sra2bam.sh

## samtools quickcheck to validate bams

# list of bams to be filtered 
ls *.bam >> bamlist.txt
# bring bed file to current directory
cp ${homeDir}/Ne-sQTL/data/12-07-2018/GRCh37.bed $PWD
# filter unplaced contigs
jid2=$(sbatch --wait --dependency=afterok:${jid1##* } ${homeDir}/Ne-sQTL/src/12-14-2018/filter_bam.sh)

## samtools quickcheck

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

## make it v7 Analysis freeze juncs ONLY && GENOTYPED

#ls SRR* >> juncfiles.txt
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



## Step 5 - QTLtools Preparation
################################################


# prepare files for QTLtools
ls *qqnorm*.gz >> leafcutterphenotypes.txt 
# important: render these files compatible with QTLtools
sbatch --wait ${homeDir}/Ne-sQTL/src/12-18-2018/QTLtools-Filter.sh
ls *.qtltools >> qtltools-input.txt
# generate the corresponding tbi files
rm Ne*tbi
for i in {1..22}; do tabix -p bed Ne-sQTL_perind.counts.gz.qqnorm_chr${i}.gz.qtltools; echo "Bedding chromosome $i"; done

mv ${homeDir}/Ne-sQTL/data/01-22-2019/GTExTissueKey.csv $PWD
# get the tissue sites for each corresonding sra file
Rscript --vanilla ${homeDir}/Ne-sQTL/src/01-09-2019/sraTissueExtract.R ${homeDir}/Ne-sQTL/data/SraRunTable.txt GTExTissueKey.csv


# submit each LF phenotype file to sraNameChangeSort as command line variable as well as tissue_table.txt
for phen in *qqnorm*.gz.qtltools; do Rscript ${homeDir}/Ne-sQTL/src/01-15-2019/sraNameChangeSort.R $phen tissue_table.txt ; done
rm *Leukemia*

#this code is problematic - I can't have any whitespaces in filenames, yet the that's all they have in tissue_table
# getting the tissue names from metadata - 48 tissues in all
cat tissue_table.txt | cut -f3 | awk '{if(NR>1)print}' |  awk '!seen[$0]++' > tissuenames.txt
#mkdir cattissues
mkdir WHLBLD/


# concatenate all tissues from chr 1-22 such that all chromosomes are found in one file per tissue
#for i in {1..48}; do line=`sed "${i}q;d" tissuenames.txt`; echo "Concatenating $line..."; for q in {1..22}; do echo "Chr $q..."; awk 'FNR==1 && NR!=1{next;}{print}' "${q}_${line}.txt" > "${line}".txt; rm "${q}_${line}.txt" done; done #formerly cat "${q}_${line}.txt" >> ... let's see if this works ; did this to prevent writing the header multiple times


#specifically just for whole blood samples
# compile phenotype file from chromosome-specific files
# print the header once

ml htslib
head -1 1_WHLBLD.txt > phen_fastqtl.bed
# print each chromosome file without the header
for i in {1..22}
do
  cat ${i}_WHLBLD.txt | sed -e1,1d >> phen_fastqtl.bed
done

ml bedtools
# INSERT CODE TO SWAP SRR FOR GTEX-SUBJECT ID HERE
bedtools sort -header -i phen_fastqtl.bed > WHLBLD.pheno.bed
bgzip WHLBLD.pheno.bed
tabix -p bed WHLBLD.pheno.bed.gz
rm phen_fastqtl.bed

mv *_WHLBLD.txt WHLBLD/

# Parts of this next line may appear redundant. That's okay. it just creates a directory for each concatenated tissue phenotype file and moves it in there.
#for file in cattissues/*; do i=$(echo $file); q=$(basename "$i"); filename="${q%.*}"; filename=$(echo ${filename#*_}); mkdir tissues/"${filename}"; mv "$file" tissues/"${filename}"; done
# at this point, you want to pass each tissue PC file and the leafcutter PC file as command-line arguments into an R script that concatenates the PCs by GTEX ID

# download genotype covariates
wget https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL_covariates.tar.gz
tar -xzf GTEx_Analysis_v7_eQTL_covariates.tar.gz


Rscript --vanilla ${homeDir}/src/01-02-2019/mergePCs.R ../Ne-sQTL_perind.counts.gz.PCs Whole_Blood.v7.covariates.txt ../tissue_table.txt

echo "Concatenating Whole Blood covariates..."

# this code is an absolute mess. I need to clean it up.

## Step 4 - Mapping sQTLs using QTLtools
################################################


##### Make this part useable for any tissue and not just whole blood


#for loop for QTLtools nominals
for i in {1..100}; do ./QTLtools_1.1_Ubuntu14.04_x86_64 cis \
	--vcf  GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz \
	--bed "WHLBLD.pheno.bed.gz" \
	--cov  "Whole_Blood.v7.covariates_output.txt" \
	--nominal 1  \
	--chunk $i 100 \
	--out "WHLBLD_nominals_chunk_${i}.txt"
done

cat WHLBLD_nominals_chunk_*.txt | gzip -c > nominals.all.chunks.txt.gz
#sbatch --export=tissue=WHLBLD.txt.gz QTLtools-NomPass.sh

ls WHLBLD_* | sort -V >> WHLBLD_chunks.txt

sbatch --wait NomPassExtractCall.sh

cat WHLBLD_nominals_chunk_*_out.txt | gzip -c > nominals.all.chunks.NE_only.txt.gz

sbatch --wait PermPass.sh

cat WHLBLD_permutations_chunk_*.txt | gzip -c > permutations_full.txt.gz

Rscript ../../../../../progs/QTLtools/script/runFDR_cis.R permutations_full.txt.gz 0.05 permuatations_full_FDR

cat WHLBLD_conditional_chunk_*.txt | gzip -c > conditional_full.txt.gz

echo "done"
