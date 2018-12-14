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

homeDir = pwd -P
# Step 1 - Convert .sra to .bam files
cd ncbi/sra
# store all .sra names into text file for job array
ls *.sra >> sralist.txt
# submit batch job
sbatch sra2bam.sh