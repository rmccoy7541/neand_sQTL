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
* LeafCutter 0.2.8
* QTLtools
* SAMtools or HTSlib
* SRA-Toolkit
* GTEx RNA-Seq SRAs downloaded using SRA-Toolkit's `prefetch`
* GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz
 downloaded using SRA-Toolkit's `prefetch` and Aspera
	* SHA1SUM: 723684f4bc6dac3c7d14e0aa27069be7f958d5b7
	* Downloaded Tue 11 Dec 2018 6:44:34 PM EST
* GTEx WGS VCF TBI (index) generated using `htslib`'s `tabix -p vcf`
	* SHA1SUM:  
* GTEx SRA Run Table from dbGaP
	* SHA1SUM: 5482b35ec85d10a9a378a0d4915f61a783a6d09a
	* Downloaded Sun 30 Dec 2018 05∶29∶53 PM EST


### Installing
How to install what they will need to install.

## Running the Tests
Step-by-step instructions with intermediate ends and full explanations of important code.

## Built With
* R 3.4.4 (2018-03-15) -- "Someone to Lean On"
* Python 2.7-anaconda
* Bash
* MARCC - Maryland Advanced Research Computing Center
* LeafCutter
* 

## Contributing
* Thank the important people
* The MARCC guys

## Authors

## Acknowledgements
	
## Bibliography

1. Rogers Ackermann, Rebecca & Mackay, Alex & L. Arnold, Michael. (2015). The Hybrid Origin of “Modern” Humans. Evolutionary Biology. 43. 10.1007/s11692-015-9348-1. 
2. Silver, L. (2001). QTL (Quantitative Trait Locus). Encyclopedia of Genetics, 1593-1595. doi:10.1006/rwgn.2001.1054
