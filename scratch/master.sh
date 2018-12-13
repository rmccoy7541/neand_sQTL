#!/bin/bash

# This is the master script that prepares and submits all jobs for LeafCutter
# This script uses sra-toolkit binary files, so at the highest directory, it is assumed that the files are formatted as such:
# sra-toolkit/		ncbi/		leafcutter/		<CartFile>.krt		NE-sQTL/		master.sh
# It is also assumed that the GTEx VCF file for the whole genome has already been downloaded (see Documentation for details)
# Finally, please make sure that you have also downloaded the SRA files (again, please see Documentation for details)

homeDir = $PWD
# Step 1 - Convert .sra to .bam files
cp sra-toolkit/bin/sam-dump ncbi/sra
cd ncbi/sra
# store all .sra names into text file for job array
ls *.sra >> sralist.txt
# submit batch job
