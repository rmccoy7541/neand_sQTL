#!/bin/bash

# This is the master script that prepares and submits all jobs for LeafCutter
# At the highest directory, it is assumed that the files are arranged as such:
# ncbi/		leafcutter/		NE-sQTL/		master.sh
# It is also assumed that the GTEx VCF file for the whole genome has already been downloaded and resides in ncbi/files/ 
# 	(see Documentation for details)
# Finally, please make sure that you have also downloaded the SRA files in ncbi/sra
# 	(again, please see Documentation for details)

# load modules
ml samtools
ml sra-tools
ml python/2.7-anaconda
ml bcftools
ml htslib
ml
# the directory of master.sh
homeDir = $(pwd -P)
# convert .sra to .bam files
cd ncbi/sra
# store all .sra names into text file for job array
ls *.sra >> sralist.txt
# submit batch job, return stdout in $RES
jid1=$(sbatch --wait ${homeDir}/NE-sQTL/src/12-14-2018/sra2bam.sh)
# list of bams to be filtered 
ls *.bam >> bamlist.txt
# bring bed file to current directory
cp ${homeDir}/NE-sQTL/data/12-07-2018/GRCh37.bed $PWD
# filter unplaced contigs
jid2=$(sbatch --wait --dependency=afterok:${jid1##* } ${homeDir}/NE-sQTL/src/12-14-2018/filter_bam.sh)
ls *.filt >> filtlist.txt
# maybe inclue an if-statement after each sbatch that would catch any non-zero exit codes and abort the program
sbatch --wait $(homeDir)/NE-sQTL/src/12-14-2018/rename_gtex.sh
ls GTEX* >> gtexlist.txt
mkdir juncfiles
sbatch --wait --dependency=afterok:${jid2##* } ${homeDir}/NE-sQTL/src/12-10-2018/bam2junccall.sh
mv *.junc juncfiles/
# include a command here that would strip the .junc extension from our boys
cd juncfiles/
# strip junc files
for junc in $PWD/*; do mv ${junc} ${junc%.*}; echo ${junc%.*}; done
# put all of the renamed junc files in a text
ls GTEX* >> test_juncfiles.txt
# intron clustering
mkdir intronclustering/
python $homeDir/leafcutter/clustering/leafcutter_cluster.py -j test_juncfiles.txt -r intronclustering/ -m 50 -o testNE_sQTL -l 500000
#PCA calculation
python $homeDir/leafcutter/scripts/prepare_phenotype_table.py testNE_sQTL_perind.counts.gz -p 10 # works fine; there's no option for parallelization.
# indexing and bedding
ml htslib; sh testNE_sQTL_perind.counts.gz_prepare.sh
# filter the genotype file using bcftools
sbatch --wait ${homeDir}/NE-sQTL/src/12-18-2018/bcf_tools.sh $homeDir
# index our friend with tabix
tabix -p vcf GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz
# call FastQTL
ls *qqnorm*.gz >> leafcutterphenotypes.txt 
sbatch --wait ${homeDir}/NE-sQTL/src/12-18-2018/FastQTL.sh



