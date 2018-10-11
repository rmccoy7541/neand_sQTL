# neanderthal-sqtl

## Introduction
This is my (Arta Seyedian) first project with Dr. Rajiv McCoy's group. Today is 9/26/2018 and I am applying principles of documentation and organziation to this project repository. Most of the work I have done so far 
has been without this method of organization so please forgive me if it's a bit hard to go through these directories and files and fully understand what I have been going for here.

The crux of this project is to find sQTLs (splicing quantitative trait loci) among the DNA modern humans have inherited from neanderthals as a result of interbreeding that took place between 100,000 - 40,000 years ago [1]. 
These sQTLs basically determine mRNA transcript splicing patterns, which lead to different proportions of mRNA isoforms which leads to differential protein expression, which impacts the expression of quantitative traits,
such as height [2].
========================================================================================================================================
### As of 10/11/2018
The repo was cloned directly onto MARCC and reorganized the data files and purged some useless stuff, including some pre-concatenated output files. Best thing I did was find a bash command that finds and appends to the gitignore file all files greater than 100MB, so I can now host this project directly on MARCC and update it through MARCC without overcomplicating getting stuck accidentally committing a massive file.

### As of 10/09/2018
I am writing this on 10/11 because I forgot to update this lab journal with what I did. I finished writing the R script that would merge the altrans tissue sQTL data with the Neanderthal rsID data and also the Bash script that would call that R script. Generated merged files for each of the Neanderthal regions and each of the tissue types. Stored them in files; now they're in /data. 

### As of 10/04/2018
I wrote some code in /scratch that merges one altrans tissue sQTL file with one neanderthal rsID file, while stripping some of the unnecesary columns. The next step is to figure out how to automate this for each one of the 3 neanderthal files, such that one neanderthal file (let's say, EASmatch.txt) merges with each one of the dozen or so tissue data files and either appends this merge to the resulting data table from the first merge, and includes maybe some header information about which tissue file that it merged with OR (and this is probably what I'll end up doing because it just occurred to me and makes more sense and seems obvious in hindsight) just make a merged file for each one of the tissue files that the neanderthal rsIDs merged with.

### As of 9/26/2018
I'm still at the beginning stages of the project. No serious analysis has taken place just yet; however, I have matched the neandertal SNP positions with the human rsIDs from the 1k data. The next step would be to match
these neanderthal rsID data with known sQTLs sequenced by altrans as a proof of concept, before we move onto using LeafCutter on GTEx data to create a full map of human sQTL data to match against the neanderthal rsIDs.

========================================================================================================================================
	
## Bibliography

1. Rogers Ackermann, Rebecca & Mackay, Alex & L. Arnold, Michael. (2015). The Hybrid Origin of “Modern” Humans. Evolutionary Biology. 43. 10.1007/s11692-015-9348-1. 
2. Silver, L. (2001). QTL (Quantitative Trait Locus). Encyclopedia of Genetics, 1593-1595. doi:10.1006/rwgn.2001.1054
