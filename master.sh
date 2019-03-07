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
homeDir=$(echo ~/work/)
scripts=$(echo /home-1/aseyedi2@jhu.edu/work/aseyedi2/neand_sQTL/src/primary/)
data=$(echo /home-1/aseyedi2@jhu.edu/work/aseyedi2/neand_sQTL/data/)
ncbiFiles=$(echo /scratch/groups/rmccoy22/Ne_sQTL/files/)

## Step 1a - Conversion & Validation
################################################
# convert .sra to .bam files
cd Ne_sQTL/sra
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
# No longer renaming SRAs until after leafcutter
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


## Step 4 - Genotype & Covariates Preparation
################################################
# filter the non-biallelic sites from genotype file using bcftools; see script for details
sbatch --wait ${scripts}/sh/bcf_tools.sh $homeDir
# index our friend with tabix
echo "Indexing our friend..."
tabix -p vcf GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz



## Step 5 - QTLtools Preparation
################################################


# prepare files for QTLtools
ls *qqnorm*.gz >> leafcutterphenotypes.txt 
# important: render these files compatible with QTLtools
sbatch --wait ${homeDir}/sh/QTLtools-Filter.sh
ls *.qtltools >> qtltools-input.txt
# generate the corresponding tbi files
rm Ne*tbi
for i in {1..22}; do tabix -p bed Ne-sQTL_perind.counts.gz.qqnorm_chr${i}.gz.qtltools; echo "Bedding chromosome $i"; done

cp ${data}/01-22-2019/GTExTissueKey.csv $PWD
# get the tissue sites for each corresonding sra file
Rscript ${scripts}/R/sraTissueExtract.R ${data}/SraRunTable.txt GTExTissueKey.csv

# submit each LF phenotype file to sraNameChangeSort as command line variable as well as tissue_table.txt
for phen in *qqnorm*.gz.qtltools; do Rscript ${scripts}/R/sraNameChangeSort.R $phen tissue_table.txt ; done

cat tissue_table.txt | cut -f3 | awk '{if(NR>1)print}' |  awk '!seen[$0]++' > tissuenames.txt

#mkdir WHLBLD/

# figure out how to make this part generalizable, but for now make a directory for each type of tissue processed

## figure out how to make this part generalizable too
ml htslib
head -1 1_TESTIS.txt > phen_fastqtl.bed
# print each chromosome file without the header
for i in {1..22}
do
  cat ${i}_TESTIS.txt | sed -e1,1d >> phen_fastqtl.bed
done

ml bedtools

bedtools sort -header -i phen_fastqtl.bed > TESTIS.pheno.bed

bgzip TESTIS.pheno.bed

pheno=$(echo $PWD/TESTIS.pheno.bed.gz)

tabix -p bed TESTIS.pheno.bed.gz
rm phen_fastqtl.bed

mkdir sepfiles/

mv *_TESTIS.txt sepfiles/

# download genotype covariates
wget https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL_covariates.tar.gz
cov=$(tar -xzf GTEx_Analysis_v7_eQTL_covariates.tar.gz)


# Make generalizeable
Rscript ${scripts}/R/mergePCs.R Ne-sQTL_perind.counts.gz.PCs Testis.v7.covariates.txt tissue_table.txt

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

ls WHLBLD_* | sort -V >> WHLBLD_chunks.txt

#Extract Neanderthal sequences
cp ${data}02-11-2019/tag_snps.neand.EUR.bed $PWD
sbatch --wait --export=listPath=$PWD,tissue=$(echo TESTIS),scripts=$scripts ${scripts}sh/NomPassExtractCall.sh

cat WHLBLD_nominals_chunk_*_out.txt | gzip -c > nominals.all.chunks.NE_only.txt.gz

mkdir nominals; mv *_nominals_* nominals/

#Call permuatation pass
sbatch --wait --export=VCF=$ncbiFiles/GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz,pheno=$pheno ${scripts}/sh/PermPass.sh

cat WHLBLD_permutations_chunk_*.txt | gzip -c > permutations_full.txt.gz

mkdir permutations; mv *_permutations_* permutations/

Rscript ${homeDir}/progs/QTLtools/script/runFDR_cis.R permutations_full.txt.gz 0.05 permuatations_full_FDR

sbatch --wait --export=VCF=$ncbiFiles/GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz,pheno=$pheno ${scripts}/sh/CondPass.sh

cat WHLBLD_conditionals_chunk_*.txt | gzip -c > conditional_full.txt.gz

echo "done"
