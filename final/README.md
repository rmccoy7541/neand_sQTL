# Neandertal-sQTL
### Finding Neandertal-introgressed splice-altering variation

### Introduction
Genetic variation influencing patterns of pre-mRNA splicing is increasingly recognized as a key source of phenotypic variation and disease risk among contemporary humans. Little research has been done, however, into the role of splicing variation in hominin evolution. In our study, we quantified splice-altering effects of introgressed archaic alleles in human gene expression data by employing the use of LeafCutter, SPrime, and QTLTools.

### Getting Started
1. This guide assumes that the user has already downloaded the .CRAM and .CRAI GTEx v7 analysis freeze sample files as well as the GTEx v7 genotype file.

2. You must also first filter the VCF of all non-biallelic sites. See `filter_vcf.sh` for details. 

3. Additionally, you must make a modification to the LeafCutter's `scripts/prepare_phenotype_table.py`, line 84:
```
        # If ratio is missing for over 40% of the samples, skip
        if tmpvalRow.count("NA") > len(tmpvalRow)*0.4:
            continue
```

	By default, LeafCutter looks for intron clusters with excision ratios present in over 60% of the samples. However, we want the broadest possible scope, so we include all intron clusters regardless of how prevelant they are throughout the samples.


```
        # If ratio is missing for over 100% of the samples, skip
        if tmpvalRow.count("NA") > len(tmpvalRow)*1:
            continue
```

4. **Lastly, set the directory variables for your project. They can be found at the top of the master script.** This is vital and failure to do so will prevent the pipeline from working.

#### Prerequisites
* LeafCutter 0.2.8
* QTLtools (v1.1)
* HTSlib (v1.9)
* BCFtools (v1.9)
* SRA-Toolkit (v2.9.6-1)
* Genome Analysis Toolkit (v4.1.3.0)
* GTEx WGS VCF downloaded using SRA-Toolkit's prefetch and Aspera
	* SHA1SUM of phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1.GRU.tar: 3ef252d20b52da537de95bc00ca42c49093771eb
	* Downloaded Tue 11 Dec 2018 6:44:34 PM EST
* High performance cluster running with SLURM (v19.05)	

## Running the Tests
Step-by-step instructions with intermediate ends and full explanations of important code.

### Filter Bi-allelic Sites from 

### Converting CRAM to JUNC
This step we performed on an Amazon Web Services server. `src/AWS/make_junc.sh` will convert your CRAM files to JUNC so long as you provide the correct directory to `leafCutterDirectory` in `bam2junc.sh`.

### JUNC Intron Clustering
Once you have all of your JUNC files, transfer them from AWS to your SLURM system and into a working directory. Obtain the path of your working directory and use it to set `$basewd` in `master.sh`. Call `master.sh` **after you make sure that the project variables at the top of the file are filled out with the correct paths.**

## Built With
* R 3.4.4 (2018-03-15) -- "Someone to Lean On"
* Python 2.7-anaconda
* Bash
* MARCC - Maryland Advanced Research Computing Center
* LeafCutter 0.2.8


## Contributing
* Thank the important people
* The MARCC guys

## Authors

## Acknowledgements
