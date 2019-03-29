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
scripts=$(echo /home-1/aseyedi2@jhu.edu/work/aseyedi2/neand_sQTL/src/primary/)
# data dir
data=$(echo /home-1/aseyedi2@jhu.edu/work/aseyedi2/neand_sQTL/data/)
# ncbi/files/
ncbiFiles=$(echo /scratch/groups/rmccoy22/Ne_sQTL/files/)
# IF YOU ALREADY HAVE NON-BIALLELIC INDEXED VCF
VCF=$(echo /scratch/groups/rmccoy22/Ne_sQTL/files/GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz)
# input directory with sra files here
sra=$(/home-1/aseyedi2@jhu.edu/work/Ne_sQTL/sra/lung_skinEx_thy)
# leafcutter directory here
leafCutter=$(echo /scratch/groups/rmccoy22/aseyedi2/leafcutter)

## Step 1a - Conversion & Validation
################################################
# convert .sra to .bam files
cd $sra
# store all .sra names into text file for job array
ls *.sra >> sralist.txt
# sra2bam, most computationally intensive step
sbatch --wait --export=sraListPath=$PWD,homeDir=$homeDir ${scripts}/sh/sra2bam.sh

## samtools error check, remove broken bams
samtools quickcheck *bam 2> samtools_err_bam.txt

cat samtools_err_bam.txt | cut -d'.' -f1,2,3 > failedbams.txt

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
echo "Filtering unplaced contigs..."
sbatch --wait ${scripts}/sh/filter_bam.sh

ls *.filt >> filtlist.txt

samtools quickcheck *filt 2> samtools_err_filt.txt

cat samtools_err_filt.txt | cut -d'.' -f1,2,3 > failedfilt.txt

for i in $(cat failedfilt.txt)
do
   echo "$i is broken, removing now..." >> log
   rm $i
done


## Step 2 - Intron Clustering
################################################
mkdir juncfiles
## IMPORTANT: changed leafcutter's bam2junc.sh to directly call
## the bin for samtools
echo "Converting bam files to junc..."
sbatch --wait ${scripts}/sh/bam2junccall.sh
mv *.junc juncfiles/
cd juncfiles/
# strip junc files - STILL WITH RUN ID 'SRR######'
find -type f -name '*.sra.bam.filt.junc' | while read f; do mv "$f" "${f%.sra.bam.filt.junc}"; done

# put all of the renamed junc files in a text
ls SRR* >> juncfiles.txt
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
for i in 1_*.txt; do echo $i | cut -d'_' -f 2| cut -d'.' -f 1 > tissuesused.txt; done

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
for line in $(cat GTExCovKey.csv)
do
   full=$(echo $line | awk -F',' '{print $1}')
   abb=$(echo $line | awk -F',' '{print $2}')
   if grep "$abb" tissuesused.txt; then
      # Nominal Pass
      sbatch --wait --export=VCF=$VCF,pheno=$(echo $abb/$abb.pheno.bed.gz),tissue=$(echo $abb/$abb),covariates=$(echo $abb/$full.v7.covariates_output.txt) ${scripts}/sh/NomPass.sh
      
      for i in {1..100}; do

         cat $abb/${abb}_nominals_chunk_${i}.txt | gzip -c >> $abb/$abb.nominals.all.chunks.txt.gz

      done

      ls $abb/$abb_* | sort -V >> $abb/${abb}_chunks_list.txt

      #Extract Neanderthal sequences
      cp ${data}02-11-2019/tag_snps.neand.EUR.bed $abb
      sbatch --wait --export=listPath=$PWD/$abb,tissue=$(echo $abb),scripts=$scripts ${scripts}sh/NomPassExtractCall.sh
      rm $abb/tag_snps.neand.EUR.bed

      for i in {1..100}; do
         cat ${abb}/${abb}_nominals_chunk_$i_out.txt | gzip -c >> ${abb}/$abb.nominals.all.chunks.NE_only.txt.gz
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

      mkdir ${leafCutter}/../sQTL/$abb

      Rscript ${scripts}/R/QQPlot-Viz.R /home-1/aseyedi2@jhu.edu/work/aseyedi2/sQTL/$abb $abb/$abb.nominals.all.chunks.NE_only.txt.gz $abb/$abb.permutations_full.txt.gz ${data}02-11-2019/tag_snps.neand.EUR.bed
   fi
done
echo "Done"
exit
