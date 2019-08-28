# Neandertal-sQTL
### Finding Neandertal-introgressed splice-altering variation

### Introduction
Genetic variation influencing patterns of pre-mRNA splicing is increasingly recognized as a key source of phenotypic variation and disease risk among contemporary humans. Little research has been done, however, into the role of splicing variation in hominin evolution. 

### Getting Started
1. Since GTEx v8 includes a splicing QTL analysis, the pipeline only includes the downloading and analysis of the splicing files provided by GTEx. In order to do this, you must have a valid dbGaP login and Google Cloud account with the ability to download the necessary files (`gtex_resources/GTEx_Analysis_v8_sQTL_all_associations/`). This will require a large amount of storage, which we are assuming you have on your cluster. In the case of Google Cloud, you will be billed for your usage.

2. Go to dbGaP and download `phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1.GRU.tar` by generating the `.krt` file and using `sra-toolkit`'s `prefetch`. Here's how we did it:
```
prefetch --max-size 1000G < cart file > ; tar -xvf < kart >
```

3. Place this file in the project directory by first cloning this project onto your desired partition:
```
git clone https://github.com/mccoy-lab/neandertal-sQTL.git
```
Then, simply move the untarred VCF subdirectory (`phg001219.v1.GTEx_v8_WGS.genotype-calls-vcf.c1.GRU/`) into the top-level of the project for execution.

#### Prerequisites
* SRA-Toolkit (v2.9.6-1)
* High performance cluster running with SLURM (v19.05)	

## Instructions

## Built With
* R 3.4.4 (2018-03-15) -- "Someone to Lean On"
* Bash
* MARCC - Maryland Advanced Research Computing Center

## Authors
- Dr. Rajiv McCoy, PhD; Assistant Professor of Biology at Johns Hopkins University
	- https://mccoy-lab.org/
- Arta Seyedian, Research Technologist at Johns Hopkins University

## Acknowledgements
* Thank you to the kind folks over at Maryland Advanced Research Computing Center (MARCC) for helping us with the logistics of computing these results.
* Thank you to the other members of the McCoy lab for their intellectual and moral support.
	* Joel Espinoza
	* Katie Farney
	* Natalie Murphy
	* Margaret Starostik
	* Stephanie Yan
