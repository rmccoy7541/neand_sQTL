#!/bin/bash

# This is the master script that prepares and submits all jobs for LeafCutter
# At the highest directory, it is assumed that the files are formatted as such:
# ncbi/		leafcutter/		<CartFile>.krt		NE-sQTL/		master.sh
# It is also assumed that the GTEx VCF file for the whole genome has already been downloaded (see Documentation for details)
# Finally, please make sure that you have also downloaded the SRA files (again, please see Documentation for details)

# load modules
ml samtools
ml sra-tools
ml python/2.7-anaconda
ml

homeDir = $(pwd -P)
# Step 1 - Convert .sra to .bam files
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
sbatch --wait $(homeDir)/NE-sQTL/src/12-14-2018/rename_gtex.sh
ls GTEX* >> gtexlist.txt
mkdir juncfiles
sbatch --wait ${homeDir}/NE-sQTL/src/12-10-2018/bam2junccall.sh
mv *.junc juncfiles/
# include a command here that would strip the .junc extension from our boys
