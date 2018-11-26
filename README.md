# neanderthal-sqtl

## Introduction
	This is my (Arta Seyedian) first project with Dr. Rajiv McCoy's group. Today is 9/26/2018 and I am applying the general guidelines for documentation and organziation to this project repository. Most of the work I have done so far 
has been without this method of organization so please forgive me if it's a bit hard to go through these directories and files and fully understand what I have been going for here.

	The crux of this project is to find sQTLs (splicing quantitative trait loci) among the DNA modern humans have inherited from neanderthals as a result of interbreeding that took place between 100,000 - 40,000 years ago [1]. 
These sQTLs basically determine mRNA transcript splicing patterns, which lead to different proportions of mRNA isoforms which leads to differential protein expression, which impacts the expression of quantitative traits,
such as height [2].


## Getting Started
	Placeholder for ultimately instructing users on how to run all of these scripts to get identical output.

### Prerequisites
	* Access to a high-performance cluster
	* LeafCutter 0.2.8
	* FastQTL v2.184
	* SAMtools or HTSlib

### Installing
	How to install what they will need to install.

## Running the Tests
	Step-by-step instructions with intermediate ends and full explanations of important code.

## Built With
	* R 3.4.4 (2018-03-15) -- "Someone to Lean On"
	* Python 2.7-anaconda
	* Obviously some other stuff

## Contributing
	* Thank the important people

## Authors

## Acknowledgements


Updates
----------------------------------------------------------------------------------------------------------------------------------------
### 11/26/2018
	The more things change, the more they stay the same. We are meeting with the director of MARCC on Wednesday to discuss our options when dealing with big ol' datasets. Looks like using `~/scratch` is a no-go, but `~/work` apparently has up to 50TB of space AND is a scratch partition, so that's what we're going to use. We're still waiting to meet with the dude on Wednesday but in the meantime, I'm going to play around with a small subset of the total data to see what's what. I already made some progress last week; I came up with a little ditty like this: `for f in $PWD/*.sra; do ./sam-dump $f | samtools view -bS > $f.bam; done`. However, this one-off command doesn't utilize parallelization, so Rajiv recommended something that makes use of SLURM's capabilities:
```#!/bin/bash
#SBATCH --job-name=convert_gtex
#SBATCH --time=12:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-10%10
#SBATCH --output=/scratch/users/rmccoy22@jhu.edu/logs/slurm-%A_%a.out

#### load and unload modules you may need
module load samtools

# test.txt contains a list of dbGaP run IDs (e.g., SRR1068687) 
line=`sed "${SLURM_ARRAY_TASK_ID}q;d" /scratch/users/rmccoy22@jhu.edu/test.txt`

echo "Text read from file: $line"

cd /scratch/users/rmccoy22@jhu.edu/dbgap/sra

# Download the .sra file according to SRA ID
~/progs/sratoolkit.2.9.2-centos_linux64/bin/prefetch \
  --ascp-path '/software/apps/aspera/3.7.2.354/bin/ascp|/software/apps/aspera/3.7.2.354/etc/asperaweb_id_dsa.openssh' $line

# Convert from .sra to .bam
~/progs/sratoolkit.2.9.2-centos_linux64/bin/sam-dump $line.sra |\
  samtools view -bS - > $line.bam
#mv $line.bam /work-zfs/rmccoy22/rmccoy22/gtex/RootStudyConsentSet_phs000424.GTEx.v7.p2.c1.GRU/BamFiles/$line.bam
#rm /scratch/users/rmccoy22@jhu.edu/dbgap/sra/$line.*

echo "Finished with job $SLURM_ARRAY_TASK_ID"```

<<<<<<< HEAD
	Since we won't need to delete and move things back and forth, I have commented-out the parts of the code that are now both practically and conceptually obsolete. I'm going to use the following command: `nohup ./prefetch -X 50000000000 /scratch/groups/rmccoy22/sratools/cart_prj19186_201811221749.krt &`
=======
	Since we won't need to delete and move things back and forth, I have commented-out the parts of the code that are now both practically and conceptually obsolete.
>>>>>>> master

	I am using the following command to download a subset of data for use: `nohup ./prefetch -X 50000000000 /scratch/groups/rmccoy22/sratools/cart_prj19186_201811221749.krt &`

### 11/12/2018
	Our journey continues. Rajiv accidentally deleted the dbGaP folder so now we have to download the sra files again. The cart file I was using before was incomplete so I got a new one. What I need to do now is to download the sra files in chunks to /scratch, convert the files to .bam, then move the converted files to /data. A good thing to keep in mind is to set the max size for prefetch to a high number so that the download won't be capped. I used the command `./prefetch -s ~/scratch/ncbi/cart_prj19186_201811121649.krt > ../../sizeOfCart.txt` to get a text file with the sizes of all the sra files so that I can find the best possible chunking. Actually, now that I think about it, that's not necessary; chunking by 3TB is probably fine and that's what I'm going to do, but it's cool that I have the sizes of the sra files. So now I have to write a bash script that will download to scratch in a 3TB chunk, convert to *.bam, move *.bam to /data, delete pre-converted data as well as bam data on /scratch, repeat for remainder of sra files in 3TB chunks.

### 11/09/2018
	Disregard that last entry, the download was unsuccessful and after speaking with Rajiv, he found a repository of GTEx data he had downloaded previously minus some larger files had yet to download. He offered to convert the files to `.bam` format for me after repeated difficulties on my end working with sra-toolkit. The total size of the RNA-Seq data that will be used is ~16TB. Additonally, he will supply me the required `.vcf` file for QTL mapping via FastQTL. I did some investigation into potentially creating an R shiny app to visualize the data and it seems not only easy, but clear cut for LeafCutter users.

### 11/08/2018
	After weeks of learning how to use leafcutter, installing a dependent QTL mapper (FastQTL), and fixing certain problems with using the GTEx data, I was finally approved by Rajiv as a downloader on dbGaP, and am about to download RNA-Seq GTEx data from dbGaP using sra-toolkit's `prefetch` command. The command will look something like this: `nohup ./prefetch -v -O /home-1/aseyedi2@jhu.edu/data/aseyedi2/GTEx_SRA cart_DAR72305_201811081445.krt &` and I will be logged into the dtn2 server, which is intended for large data transfers. It's about 16 TB. 

### 11/02/2018
	Got rid of the contents of the scratch folder, updated the README for `TissueMerge.R`. Added analysis folder.

### As of 10/17/2018
	Fixed a problem with the column headers in the merged Neanderthal x Altrans data.

### As of 10/15/2018
	I realized that I didn't include the Ensembl ID codes in my Neanderthal x Altrans data merge, so changed the TissueMerge.R code to make it do that.

### As of 10/11/2018
	The repo was cloned directly onto MARCC and reorganized the data files and purged some useless stuff, including some pre-concatenated output files. Best thing I did was find a bash command that finds and appends to the gitignore file all files greater than 100MB, so I can now host this project directly on MARCC and update it through MARCC without overcomplicating getting stuck accidentally committing a massive file.

### As of 10/09/2018
	I am writing this on 10/11 because I forgot to update this lab journal with what I did. I finished writing the R script that would merge the altrans tissue sQTL data with the Neanderthal rsID data and also the Bash script that would call that R script. Generated merged files for each of the Neanderthal regions and each of the tissue types. Stored them in files; now they're in /data. 

### As of 10/04/2018
	I wrote some code in /scratch that merges one altrans tissue sQTL file with one neanderthal rsID file, while stripping some of the unnecesary columns. The next step is to figure out how to automate this for each one of the 3 neanderthal files, such that one neanderthal file (let's say, EASmatch.txt) merges with each one of the dozen or so tissue data files and either appends this merge to the resulting data table from the first merge, and includes maybe some header information about which tissue file that it merged with OR (and this is probably what I'll end up doing because it just occurred to me and makes more sense and seems obvious in hindsight) just make a merged file for each one of the tissue files that the neanderthal rsIDs merged with.

### As of 9/26/2018
	I'm still at the beginning stages of the project. No serious analysis has taken place just yet; however, I have matched the neandertal SNP positions with the human rsIDs from the 1k data. The next step would be to match
these neanderthal rsID data with known sQTLs sequenced by altrans as a proof of concept, before we move onto using LeafCutter on GTEx data to create a full map of human sQTL data to match against the neanderthal rsIDs.

	
## Bibliography

1. Rogers Ackermann, Rebecca & Mackay, Alex & L. Arnold, Michael. (2015). The Hybrid Origin of “Modern” Humans. Evolutionary Biology. 43. 10.1007/s11692-015-9348-1. 
2. Silver, L. (2001). QTL (Quantitative Trait Locus). Encyclopedia of Genetics, 1593-1595. doi:10.1006/rwgn.2001.1054
