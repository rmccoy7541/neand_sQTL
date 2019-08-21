#!/bin/bash
# Snakefile wrapper
# modify this to load the modules most appropriate to your scheduler as necessary
module load bedtools
module load samtools
module load sra-tools
module load python/2.7-anaconda
module load bcftools
module load htslib
module load R
module load gcc

#TODO figure out how to call the snakefile in a way that allows for a month-long pipeline
snakemake
# https://bioinformatics.stackexchange.com/questions/4977/running-snakemake-on-cluster