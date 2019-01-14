# Updates

These updates are read from most recent date at the top to initial entry at the bottom.

NEW PIPELINE
- Generate LeafCutter phenotypes without renaming the SRAs
- Massive chr 1-22 tables with all across all samples
- Use dbGaP metadata to segregate all columns by tissue, create a new table for those samples in tissue directory
- put covariate file in tissue directories
- Change the column header in each tissue directory's folder to the full GTEX sample ID
- Run QTLtools for each tissue
also remember
- CrossRef the GTEx file that contains all of our samples of interest
- use samtools to convert cram files to bam files to use with leafcutter

### 01/14/2019
#### Mon 14 Jan 2019 10:54:43 AM EST 
~~After a couple of days of procrastinating, I am ready to face this monster. This demon. First, last night or maybe two nights ago I realized that I should probably rename these guys **before** running `prepare_phenotype.py` or whatever it's called. This is because, when it generates the covariates, it makes the Run ID the colname. I could also change the colnames of the covariate file in place, which is what I think I will do.

I'm going to try to change everything after running `prepare_phenotype.py`. This is going to be hard. I'm going to download these files so I can develop the R script locally.

Okay. I'm pretty confused and am not sure what to do. The GTEx supplied covariates are all in subject id format i.e. GTEX-.... (with an optional 5th char). I have a tissue table which has sra run IDs, corresponding tissue type, and *full* GTEx sample ID i.e. GTEX-....-....-..-.... (or something like that). Not sure how to proceed here. Do I want to change the LeafCutter PCs file to match the GTEx supplied one? My puny brain which was intentionally designed for being dumb is having a meltdown.~~

I'm confusing the hell out of myself. I think I should just roll back a step, to before `prepare_phenotype.py`.

Starting over with `python ../../../leafcutter/clustering/leafcutter_cluster.py -j juncfiles.txt -r intronclustering/ -m 50 -o NE_sQTL -l 500000`. 

Okay, I reread the notes I took. I guess I forgot that I was doing this (this being... renaming column headers and segregating based on tissue) for only the phenotype files. But I also have to do this sort of thing for the covariate files. I guess I could acheive this by simply renaming the colnames for the LC-generated PCs and then *renaming the GTEx-supplied covariate colnames* to the sample id of the tissue/subject complex.

This is so complicated. I have no idea how to make this more simple. I do not like how complicated it is. 

I'm just going to go ahead with doing this rename and rewrite thing but only for the phenotype files right now. I will later focus on the PCs.

### 01/12/2019
#### Sat 12 Jan 2019 03:29:24 PM EST 
I've done a bunch of little things and I don't feel like listing them, just check the version history. I'm now trying to set up the tissue directories. 

I feel like what I want to do right now is immensely complicated. I have a file named `tissue_table.txt` that looks like this:

```
Run	Sample_Name	body_site
SRR598484	GTEX-PW2O-0526-SM-2I3DX	Lung
SRR598124	GTEX-NPJ8-0011-R4a-SM-2HML3	Brain - Amygdala
SRR599192	GTEX-N7MT-0011-R5a-SM-2I3G6	Brain - Caudate (basal ganglia)
SRR601925	GTEX-OHPK-0526-SM-2HMJB	Lung
SRR601068	GTEX-Q2AG-0126-SM-2HMLB	Skin - Sun Exposed (Lower leg)
SRR602598	GTEX-Q2AG-0011-R9A-SM-2HMJ6	Brain - Spinal cord (cervical c-1)
SRR607586	GTEX-OXRL-0526-SM-2I3EZ	Lung
SRR608288	GTEX-OXRK-0926-SM-2HMKP	Lung
SRR600445	GTEX-Q2AG-0011-R4A-SM-2HMKA	Brain - Amygdala
SRR608344	GTEX-OIZH-0005-SM-2HMJN	Whole Blood
```
with the run id in the first column, full GTEx sample in second, and tissue in third. I also have a ton of phenotype files:

```
NE_sQTL_perind.counts.gz.qqnorm_chr10.gz.qtltools  NE_sQTL_perind.counts.gz.qqnorm_chr20.gz.qtltools
NE_sQTL_perind.counts.gz.qqnorm_chr11.gz.qtltools  NE_sQTL_perind.counts.gz.qqnorm_chr21.gz.qtltools
NE_sQTL_perind.counts.gz.qqnorm_chr12.gz.qtltools  NE_sQTL_perind.counts.gz.qqnorm_chr22.gz.qtltools
NE_sQTL_perind.counts.gz.qqnorm_chr13.gz.qtltools  NE_sQTL_perind.counts.gz.qqnorm_chr2.gz.qtltools
NE_sQTL_perind.counts.gz.qqnorm_chr14.gz.qtltools  NE_sQTL_perind.counts.gz.qqnorm_chr3.gz.qtltools
NE_sQTL_perind.counts.gz.qqnorm_chr15.gz.qtltools  NE_sQTL_perind.counts.gz.qqnorm_chr4.gz.qtltools
NE_sQTL_perind.counts.gz.qqnorm_chr16.gz.qtltools  NE_sQTL_perind.counts.gz.qqnorm_chr5.gz.qtltools
NE_sQTL_perind.counts.gz.qqnorm_chr17.gz.qtltools  NE_sQTL_perind.counts.gz.qqnorm_chr6.gz.qtltools
NE_sQTL_perind.counts.gz.qqnorm_chr18.gz.qtltools  NE_sQTL_perind.counts.gz.qqnorm_chr7.gz.qtltools
NE_sQTL_perind.counts.gz.qqnorm_chr19.gz.qtltools  NE_sQTL_perind.counts.gz.qqnorm_chr8.gz.qtltools
NE_sQTL_perind.counts.gz.qqnorm_chr1.gz.qtltools   NE_sQTL_perind.counts.gz.qqnorm_chr9.gz.qtltools
```
which each contain sample information as such:

```
[aseyedi2@jhu.edu@compute0175 intronclustering]$ zcat NE_sQTL_perind.counts.gz.qqnorm_chr17.gz.qtltools | head
#Chr	start	end	ID	.	+	SRR601068	SRR599192	SRR598484	SRR601925	SRR607586	SRR598124	SRR602598	SRR608344	SRR600445	SRR608288
17	914103	914385	17:914103:914385:clu_242_NA	.	+	0.8486631322642567	-0.9028329307690043	-0.48356192470323517	1.990047368643184	1.1073250856949786	-1.067056903921965	-0.3230811226873516	0.5399025398720513	-0.9903307055215957	-0.5089548923498618
17	914103	915086	17:914103:915086:clu_242_NA	.	+	-0.8525574183170223	0.9110166377754605	0.6842637898910064	-2.2366906721286655	-1.074273388889014	0.9684244937444385	0.39818950364007744	-0.7207118830582105	0.902606468982208	0.7055563878099839
17	915225	915928	17:915225:915928:clu_243_NA	.	+	-0.07956489889950728	0.8215433191997514	-0.7020819862644853	-1.0992940719017426	-0.06793127017856682	0.9053270737028077	0.7376160886524243	-1.9342885283098337	0.8724468436863674	0.562175939259561
17	915797	915928	17:915797:915928:clu_243_NA	.	+	0.1972239104174658	-0.7728877020815833	0.8529909148778154	1.0622766254525227	0.2988851856516097	-0.9725265742164065	-0.792124561392808	1.4717826271510717	-0.9271067287932372	-0.32466900655206166
17	1003975	1012174	17:1003975:1012174:clu_244_NA	.	+	0.8270457899025043	-1.241843521736311	0.43947833126406216	0.501251979241662	0.420143179232026	-1.1974652340294456	-0.771669548635817	1.1987005974336862	-1.1045477219149238	0.7900642238203405
17	1003975	1028518	17:1003975:1028518:clu_244_NA	.	+	-0.875979615366985	1.2027282193859805	-0.21646581467441642	-0.33516980981693184	-0.20368069440873698	1.0834411763774872	0.7201259305240991	-1.6432820168901987	1.0139816920087252	-0.6010781589979557
17	1268352	1273003	17:1268352:1273003:clu_245_NA	.	+	-1.4985649346311092	-0.22109555248029047	-0.8083123481260814	-0.00926664132977811	0.3933010555168565	0.5258359932679728	2.1669897779111342	0.07578664228527895	-0.2875599864868766	-0.12759032172880408
17	1268352	1303341	17:1268352:1303341:clu_245_NA	.	+	1.2021073155555755	0.4558006992174561	0.7615623865236025	0.44812432048529693	-0.041221266481142924	-0.4900077761500515	-2.570380323280003	-0.33421368011489794	0.47561708848328343	0.42805614716447965
17	1273035	1303341	17:1273035:1303341:clu_245_NA	.	+	-1.1873424499471859	-0.42690051379285743	-0.26782834696045893	-0.5502160519586197	0.14173083445409398	0.5220334649134223	2.1162546604968853	0.4551321272853628	-0.5747499104115467	-0.22418467313515364
```
I want to match up the column names to the tissues found in `tissue_table.txt`, make a new phenotype table for each tissue formatted the same as the original phenotype table found directly above, rename the columns to the full GTEx sample ID found in `tissue_table.txt`, and then I guess that's it. I feel like I should do this in R but it seems like a daunting task.


### 01/11/2019
#### Fri 11 Jan 2019 01:14:57 PM EST 
They're all bam files now, every single one of the, but they're still named after the SRR run number thing. I need to figure out how to rename them to the full GTEX ID e.g. `GTEX-NPJ8-0011-R4a-SM-2HML3`. ~~The problem is, they're not all uniform in the format, so I have to figure out how to detect them even though some of them are 21 characters and others 16 (without the dashes).~~ ~~Nevermind, I could use the tissue table generated by SraRunTable.txt, and match the sra to the tissue type and GTEX ID. Not sure how I would do that in bash though.~~ Actually I don't need to do any of this until after I generate the phenotype tables with leafcutter and I'm trying to excise the columns and place them in new tables in the corresponding tissue subdirectories. 

Okay I'm meeting with Rajiv about the article
- no reason to wait until project is done to start outlining different sections
	- esp intro, figures
	- flexible roadmap
- Project notebooks in R
	- good for making documentation more readable
	- include code snippets interspersed with portions of data tables
- Shiny R app
- GitHub io page step by step walkthrough
- Read that splicing review article
	- splice site qtl
		- could predict from sequence alone
	- splice enhancer/repressor sequence
		- also predict from seq alone
		- also look at sequences that control the factors themselves
- *important* it would be cool to see what splicing transcripts are completely unique to neanderthal introgression
- google docs, start entering references
- test difference between neanderthal introgressed sqtls and phenotypic variation vs non-neanderthal splice-altering mutations (or something like that)

#### Fri 11 Jan 2019 04:09:51 PM EST 
Above are the notes I took during my meeting with Rajiv about the outline of the paper. 

Anyway, on with our pipeline; newly filtered bam files are now available. 

Started conversion to junc files. Made some changes to `bam2junccall.sh`; see version history.

#### Fri 11 Jan 2019 05:31:19 PM EST 
I'm signing off for the day. I got up until the point where I need to reach into the phenotype files and match and rename headers to their full GTEx ID equivalent while sorting them into distince tissue subdirectories. **Also, I need to bear in mind that in MARCC, because of some technical difficulty that I'm having with `sra-toolkit`, the `ncbi/` file is actually `Ne-sQTL`, and in the work partition, there is no git directory. I don't know if that makes sense, but anyway, don't get confused while developing the pipeline. I doubt I'll even read this again** Good night.

### 01/09/2019
#### Wed 09 Jan 2019 02:02:37 PM EST 
Downloaded a new set of SRA files to test the pipeline with.

`tissuesite <- subset(sratabl, select=c("Run", "body_site"))`

I'm trying to figure out how to extract the tissue metadata about the samples that get downloaded. Above works, I think. I need to write it to file. ~~Not sure how I'm going to use it.~~ I'm going to use it after the leafcutter phenotype tables have been generated to figure out what each sra and/or GTEx name corresponds to. I included a little R ditty that would take care of that, `sraTissueExtract.R`. Bad name, I know.

### 01/08/2019
#### Tue 08 Jan 2019 12:00:46 PM EST 
I screwed something up. Overnight, QTLtools produced our nominal pass results, but in trying to tidy up the subdirectory names (how does "`Esophagus_Gastroesophageal_Junction.v7.covariates_output.txt_folder`" look to you?) I ended up destroying over half of the results (I used `for tissue in *; do newname=$(echo $tissue | awk -F'[_.]' '{print $1}'); mv $tissue/* $newname; done` when I *SHOULD HAVE USED* `for tissue in *; do newname=$(echo $tissue | awk -F'[.]' '{print $1}'); mv $tissue/* $newname; done`; note the missing underscore). This actually isn't that big of a deal, since I have all of the same data that I used to generate the nominal pass results, but it'll just take a long time and be a pain in the ass. That's right, I said it: "ass". There is no "undo" option in shell and it is painfully clear.

I also accidentally renamed all of the merged covariate files from `Nerve_Tibial.v7.covariates_output.txt` to just `Nerve_Tibial`, which, again, not that big of a deal. 

That does mean, however, that I no longer need to "clean things up a bit." I will have to include that lethal one liner `for tissue in *; do newname=$(echo $tissue | awk -F'[.]' '{print $1}'); mv $tissue/* $newname; done` in the master script somehow, though.

This didn't work, and according to Rajiv it's because they need to be processed by sample ID.

-covariates need sample ID
-through leafcutter, generate 22 chromosomes times however many tissues there are
-all of our analysis will be within tissues; segregate the bam files prior to using leafcutter

NEW PIPELINE
- Generate LeafCutter phenotypes without renaming the SRAs
- Massive chr 1-22 tables with all across all samples
- Use dbGaP metadata to segregate all columns by tissue, create a new table for those samples in tissue directory
- put covariate file in tissue directories
- Change the column header in each tissue directory's folder to the full GTEX sample ID
- Run QTLtools for each tissue
also remember
- CrossRef the GTEx file that contains all of our samples of interest
- use samtools to convert cram files to bam files to use with leafcutter

#### Tue 08 Jan 2019 02:56:07 PM EST 
Look at the above list. Learn to love it. We're starting over. How over? I'm downloading the samples as cram files now. Had a good meeting with Rajiv where we came to sketch out the overall structure.

At the top level we have the VCF file of the whole genome sequence. Then, we will split the analysis up into separate tissue files. We will run leafcutter on all the samples but the create a new table for each tissue that contains its corresponding phenotypes. Probably by chromosme. That sounds hard. But I'm going to do it. Probably in R. 

* VCF
	* Adipose_Subcutaneous
		* Phenotypes
		* Covariates
	* Adipose_Visceral
		* Phenotypes
		* Covariates
	* etc


I deleted everything and am starting over. Downloading new samples in cram now. UPDATE: Downloading cram files don't work. We're sticking to sra's for now.

#### Tue 08 Jan 2019 04:20:50 PM EST 
```
[aseyedi2@jhu.edu@compute0231 Ne_sQTL]$ ./prefetch -l cart_prj20712_201901081608.krt | awk -F'[||]' '{print $3}' > srafiles.txt

SRR598124
SRR598484
SRR599192
SRR600445
SRR601068
SRR601925
SRR602598
SRR607586
SRR608288
SRR608344
```
`sed -i '/^$/d' srafiles.txt`

#### Tue 08 Jan 2019 04:50:09 PM EST 
I'm redownloading the SRA files. I'll get back to this later.

### 01/07/2019
#### Mon 07 Jan 2019 02:57:28 PM EST 
**Important**: consider switching the batch scripts from `parallel` instead of `shared`, also consider using the `--exclusive` flag. 

I'm going to do some troubleshooting to see why I am being punished this way. I will do so by *okay my computer crashed* I will do so by redoing the nominal pass but naming all of the resultant files with the SLURM job metadata: `--out slurm_${SLURM_JOBID}_nominals_chr${SLURM_ARRAY_TASK_ID}_chunk${chunk}.txt` 

I also need to separate them by directory.

`for tissue in testcovmer/*; do mkdir "${tissue}_folder"; for chunk in {1..20}; do sbatch --export=tissue=$tissue,chunk=$chunk QTLtools-NomPass.sh; done; done`

```
#!/bin/bash
#SBATCH --job-name=QTLtools-NomPass
#SBATCH --time=2:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-22

# tissue is going to be $tissue and chunk $chunk
./QTLtools_1.1_Ubuntu14.04_x86_64 cis \
  --vcf GTExWholeGenomeSequenceGenotypeMatrixBiallelicOnly.vcf.gz \
  --bed testNE_sQTL_perind.counts.gz.qqnorm_chr${SLURM_ARRAY_TASK_ID}.gz.qtltools \
  --cov $tissue \
  --nominal 0.01 \
  --chunk $chunk 20 \
  --out ${tissue}_folder/slurm_${SLURM_JOBID}_nominals_chr${SLURM_ARRAY_TASK_ID}_chunk${chunk}.txt
```

also figure out how to clean things up a little bit with `$tissue` using `tissue=$(echo $tissue | awk -F'[_.]' '{print $1})'`

`# must be in testcovmer/`
`for tissue in *; do newname=$(echo $tissue | awk -F'[_.]' '{print $1}'); mv $tissue/* $newname; done`

Something like that. Maybe just figure it out later. Or not. Don't fix it if it ain't broke, right? Right.

### 01/04/2019
#### 3:09 PM 1/4/2019
Working from home today. The job array worked and now we have all 21,000+ nominal pass files ready. To merge them, I'm going to be using this shell command that I borrowed from FastQTL.
`cat *nominals_*_20.txt | gzip -c > permutations_full.txt.gz`	
...except with "nominals." There are so many files that MARCC needs a few minutes to display them.

Strangely enough, there are only a little over 15,000 output files when there should be more. I wonder why that is.


### 01/03/2019
#### Thu 03 Jan 2019 12:20:46 PM EST 
I have to figure out a way to like... tri-layer the job paralellization. I have 22 chromosomes, something like 48 different tissues with covariates, and each chromosome has 20 chunks that will be concatenated. Together, that's 22 * 48 * 20 = **21120** independent jobs that I want to run *if* I want to do this as efficiently as possible, otherwise I could just do a single job array for each chromsome that churns out nominal pass results in chunks of 1/1.

I could just make the job array only handle different chromosomes and make a nested bash for-loop that handles each tissue and chunk. For example:

```
for tissue in tissuefolder
	for chunk of 20
		do sbatch QTLtools-nompass.sh $tissue $chunk
	done
done
```
where I pass command-line variables to `QTLtools-nompass.sh`, such that, let's say the first tissue will be adipose and the first chunk 1/20, this script will call adipose chunk 1/20 for *every* chromosome 1-22 at once. Or at least I hope. Let's try it out.


`for tissue in testcovmer/*; do for chunk in {1..20}; do sbatch --export=tissue=$tissue,chunk=$chunk QTLtools-NomPass.sh; done; done`

#### Thu 03 Jan 2019 04:22:43 PM EST 
Took a while, but realized that our friend over here does **not** work.

#### Thu 03 Jan 2019 05:03:49 PM EST 
Got some advice from Rajiv. Fixed it. Was fundamentally a formatting issue. Important note: outputs all files to `testcovmer/`, which is fine I guess but not necessarily what I had in mind. Will have to deal with that in the master script.


* ~~Output the samples and sort them by tissue~~
	* ~~Get table that matches GTEX ID with tissue~~
* Run leafcutter for all samples of all tissues
* ~~Concatenate leafcutter-generated PCs with GTEx PCs~~
	* ~~Transpose and merge in R the covariates~~
* ~~Run QTLtools on each set of tissues using the concatenated PCs~~
* *Deal with the outputs going straight to testcovmer/*


### 01/02/2019
#### Wed 02 Jan 2019 01:08:58 PM EST
I'm running the QTLtools command again as a single-line type deal:
`./QTLtools_1.1_Ubuntu14.04_x86_64 cis --vcf GTExWholeGenomeSequenceGenotypeMatrixBiallelicOnly.vcf.gz --bed testNE_sQTL_perind.counts.gz.qqnorm_chr1.gz.qtltools --cov Whole_Blood.v7.covariates_output.txt --nominal 0.01 --chunk 1 1 --out WholeBlood_nominals_chr1.txt`

Right after I looked for the command, the process aborted. I keep getting this error, no matter what I do:
```
terminate called after throwing an instance of 'std::invalid_argument'
what():  stof
```
It occurs right after `Residualize phenotypes for covariates`

I'm going to talk to Rajiv about this since I tried multiple different things.

* ~~Output the samples and sort them by tissue~~
	* ~~Get table that matches GTEX ID with tissue~~
* Run leafcutter for all samples of all tissues
* ~~Concatenate leafcutter-generated PCs with GTEx PCs~~
	* ~~Transpose and merge in R the covariates~~
* Run QTLtools on each set of tissues using the concatenated PCs

#### Wed 02 Jan 2019 03:51:01 PM EST 
Okay, weird stuff: it turns out that the reason QTLtools wasn't working is because the concatenated covariate file in question (Whole blood) was missing i.e. had "NA" values for GTEX-XV7Q. I looked it up on dbGaP and our friend GTEX-XV7Q actually **has** phenotype/expression data for whole blood. So it's just not represented on the covariate file and a handful of other tissue covariate files. I emailed GTEx's help desk and in the meantime and simply going to totally exclude any samples that contain any "NA" values since that seems to be what's junking up QTLtools.

#### Wed 02 Jan 2019 05:21:26 PM EST 
Figured out how to do this good; I included removing all columns with NA's in `mergePCs.R`. Let's see if this works now.

#### Wed 02 Jan 2019 07:59:50 PM EST  
It worked. So, it's important to note: **QTLtools does not take NA values.**

### 01/01/2019
#### Tue 01 Jan 2019 04:11:17 PM EST 
Happy new year. It seems like whatever I do, I keep getting this error from QTLtools:
```
Residualize phenotypes for covariates
terminate called after throwing an instance of 'std::invalid_argument'
  what():  stof
/var/spool/slurm//job31982167/slurm_script: line 10:  3263 Aborted                 ./QTLtools_1.1_Ubuntu14.04_x86_64 cis --vcf GTExWholeGenomeSequenceGenotypeMatrixBiallelicOnly.vcf.gz --bed testNE_sQTL_perind.counts.gz.qqnorm_chr$SLURM_ARRAY_TASK_ID.gz.qtltools --cov ${tissue} --nominal 0.01 --chunk 1 1 --out ${tissue}_nominals_chr${SLURM_ARRAY_TASK_ID}.txt
```


### 12/31/2018
#### 10:39 AM 12/31/2018
It's the last day of the year and I'm trying to remember what I was working on before I took a week+ long break for the holidays. According to my chat log with Rajiv, QTLtools wasn't working when chunked, so I'm going to have to stick with only doing them in 1/1 chunks, which is fine and will take longer but it doesn't appear that I have any other choice.

Okay, MARCC isn't working with Windows right now, gonna switch computers. Just to save my train of thought, I was going to try to concatenate all of the covariate files.

#### Mon 31 Dec 2018 12:29:43 PM EST 
Concatenating the PCs worked but I ran into a problem with `QTLtools-nompass.sh`. I just realized that if I'm going to "convert" the leafcutter-generated files to a QTLtools compatible format, I'm going to need to re-index them.

Surprise, it worked. Check out the master script for details.

### 12/30/2018
#### Sun 30 Dec 2018 05:47:25 PM EST 
Christmas is over. Deal with it. I downloaded the SRA run table. It's the metadata for all of GTEx because dbGaP won't allow me to only download metadata for the runs that I want. It's in `data/`.

Okay, don't forget what you have to do:
* ~~Output the samples and sort them by tissue~~
	* ~~Get table that matches GTEX ID with tissue~~
* Run leafcutter for all samples of all tissues
* Concatenate leafcutter-generated PCs with GTEx PCs
	* Transpose and merge in R the covariates
* Run QTLtools on each set of tissues using the concatenated PCs

### 12/22/2018
#### 12:21 PM 12/22/2018

Remember to `sha1sum` all the sra's or whatever you end up downloading.

### 12/21/2018
#### Fri 21 Dec 2018 11:10:08 AM EST 
Okay, we are going to start with that script that converts all of the bed files aka the phenotype files generated by leafcutter to bedfiles compatible with QTLtools.

`zcat myFastQTLphenotypes.bed.gz | awk '{ $4=$4" . +"; print $0 }' | tr " " "\t" | bgzip -c > myQTLtoolsPhenotypes.bed.gz`

I renamed `src/12-18-2018/FastQTL.sh` to `QTLtools.sh`.

Okay, don't forget what you have to do:
* Output the samples and sort them by tissue
	* Get table that matches GTEX ID with tissue
* Run leafcutter for all samples of all tissues
* Concatenate leafcutter-generated PCs with GTEx PCs
	* Transpose and merge in R the covariates
* Run QTLtools on each set of tissues using the concatenated PCs

#### Fri 21 Dec 2018 04:57:47 PM EST 
`mergePCs.R` works.

#### Fri 21 Dec 2018 07:17:08 PM EST 
Paralellizing with QTLtools does **not** work.

### 12/20/2018
#### Thu 20 Dec 2018 01:41:29 PM EST 
* Get GTEx covariates
* ~~do PCA with QTLtools~~ I don't need to do that if I get the GTEx covariates, duh.
* Get tissue information for each SRA - either GTEx manifest or dbGAP cartinfo table
* Contact GTEx portal and ask if i want to replicate GTEx v7 eQTL analysis, which files in dbGAP would I use?
	* why are there 12k of them when that doesn't match table #?

I sent the cowboys over at GTEx an email asking them what the CRAM files were and which files they recommend for use for QTL analysis.

GTEx v7 eQTL covariates downloaded Thu 20 Dec 2018 02∶52∶59 PM EST
	739e6392c7be76adaf090bc402d23046747bcaa9  GTEx_Analysis_v7_eQTL_covariates.tar.gz


QTLtools:

Nominal pass > Permutation pass > Conditional analysis

Nominal pass tests each variant within a let's say 1 Mb range of the phenotype for an association. Then, permutation pass adjusts the p-value to account for multiple testing (all variants within range of a phenotype as well as all of the variants from all phenotype ranges; how this affects the p-value, idk). Then, conditional analysis accounts for like genetic drift and linkage disequilibrium. 

### 12/19/2018
#### Wed 19 Dec 2018 01:40:56 PM EST 
FastQTL is bunk. It doesn't work with covariates, which is why we have to use QTLtools instead. Going to read the paper on QTLtools first. Regardless of what we do, QTLtools does not accept BED files the same way FastQTL does. Here's the "quick and dirty" conversion sourced from QTLtools > Preparing Input Files: `zcat myFastQTLphenotypes.bed.gz | awk '{ $4=$4" . +"; print $0 }' | tr " " "\t" | bgzip -c > myQTLtoolsPhenotypes.bed.gz`

Gonna have to do [cis] Discover QTL in cis [nominal pass] or [cis] Discover QTL in cis [permutation pass]

### 12/18/2018
#### Tue 18 Dec 2018 10:31:49 AM EST 
	31d26b9b8db53fa92d322bf34d2a2db9a671d6f5  phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1.GRU.tar.ncbi_enc

Above is the NCBI-encrypted enormous 488 Gb VCF-TAR. I need to make a note of the shasum in the documentation. 

Okay. Our friend has been tabix'd. Let's try the FastQTL command once more.

`fastQTL -V GTExWholeGenomeSequenceGenotypeMatrixBiallelicOnly.vcf.gz -B testNE_sQTL_perind.counts.gz.qqnorm_chr1.gz -O res -L res.log --chunk 1 10`

Okay it works thankfully. Now I have to figure out how to parallelize it and also consider covariate analysis. 

`fastQTL -V GTExWholeGenomeSequenceGenotypeMatrixBiallelicOnly.vcf.gz -B testNE_sQTL_perind.counts.gz.qqnorm_chr1.gz -C testNE_sQTL_perind.counts.gz.PCs -O res -L res.log`


### 12/17/2018
#### Mon 17 Dec 2018 12:39:07 PM EST 
Okay, I've done all of the steps from making the junc files to doing PC analysis, and now I need to see if FastQTL will work. 
`fastQTL -V ../../../files/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz -B  testNE_sQTL_perind.counts.gz.qqnorm_chr1.gz -O res -L res.log --chunk 1 10`

Ok, it doesn't work, even though I renamed the bam files to their GTEX ID before generating the junc files. I'm going to talk to Rajiv.

`bcftools view -m2 -M2 -v snps --threads 23 -O z -o biallelicOnly.vcf.gz ../../../files/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz`

The above command is the filter out the non-biallelic sites. `--threads 23` because I'm running it on the dev node, and it's additional threads (dev node has 24 cores). Later, when running on MARCC, I can do something crazy like a bajillion cores probably. Afterwards, I can use tabix to split the vcf by chromosome.

#### Mon 17 Dec 2018 03:37:38 PM EST 
**I guess we would also have to figure out a way to make sure that all of the GTEX samples that were not whole genome sequenced are excluded from the FastQTL mapping somehow using bcftools. This is important since FastQTL will not work otherwise; it cannot take files with non-overlapping headers.**

#### Mon 17 Dec 2018 04:22:41 PM EST 
Remember to include a section on how to reproducibly download data. This is very important as well. Also `sha1sum`/checksum everything.

#### Mon 17 Dec 2018 05:25:51 PM EST 
It's done filtering out the non-biallelic sites. I just used `tabix -p vcf GTExWholeGenomeSequenceGenotypeMatrixBiallelicOnly.vcf.gz` to index our new friend. 

### 12/14/2018
#### 12/14/2018 10:14:38 AM
`sra2bam.sh` finished running overnight. I'm still building the master script (`master.sh`) and am figuring out how to make a pipeline by using SLURM dependencies. I have now called `filter_bam.sh` which should remove unplaced contigs. I should still remember to fire off `for file in $PWD/*.filt; do ID=$(samtools view -H $file | grep -o -m 1 'GTEX-....'); mv $file ${ID}.sra.bam.filt; done;` as a SLURM job array once this step is done. Should be easy.

#### 2:16 PM 12/14/2018
Should do something like `ls *.filt >> filtlist.txt` to make it a job array.

#### 7:18 PM 12/14/2018
Ran the junc file script after renaming each file in the entire directory. Now I have to do something like `purename="${filename%.*}"` to remove the extension from each file and rename the file accordingly, so it's just `GTEX-####`.

### 12/13/2018
#### Thu 13 Dec 2018 11:49:58 AM EST 
FastQTL isn't working. I'm going to talk to Rajiv about this.

#### Thu 13 Dec 2018 01:34:54 PM EST 
According to Rajiv, the .vcf file refers to individuals by a GTEX ID (GTEX-????). However, the corresponding phenotype files are in SRA ID format.

#### Thu 13 Dec 2018 04:59:22 PM EST  
I am the president. I am the olympics. Here is the command to fire off after filtering our SRA files: `for file in $PWD/*.filt; do ID=$(samtools view -H $file | grep -o -m 1 'GTEX-....'); mv $file ${ID}.sra.bam.filt; done;` This will find within each filtered bam file its corresponding GTEX ID and rename the file accordingly. Unfortunately, this means I will have to start the whole thing over again. Fortunately, this is a good opportunity to make the master script.

#### Thu 13 Dec 2018 05:40:41 PM EST 
Okay, so apparently `ml sra-tools` works now, so I can actually make the master script much easier now without copying the binary files into each working directory. I have submitted `sra2bam.sh` now, we'll see how it goes.

### 12/12/2018

`fastQTL -V genotypes.vcf.gz -B  phenotypes.bed.gz -O res -L res.log --chunk 1 10 `

This is the sample FastQTL line I got from the [NIH website.](https://hpc.nih.gov/apps/FastQTL.html). According to the LeafCutter, I have to use the `*qqnorm*.gz` files as the BED files. I downloaded the enormous, 447Gb genotype matrix in VCF format and now have to decode them using this command: `./vdb-decrypt -v phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1.GRU.tar.ncbi_enc GTExVCF > decryption.log`. Once that happens, I'm going to investigate fastQTL's multithreading/parallelization capabilities to see how I make a job array. 

I just learned about `wait` in Bash. I'm going to use that when writing the master script (mostly because I have to).

**I have to make a note that whomever tries to reproduce these results definitely keeps the bin files for SRAToolkit in the respective directories.**

~~Okay, the VCF is done decrypting. I need to run FastQTL on each file, and chunk the QTL mapping into 10 distinct pieces.~~

~~**BUT FIRST** I need to run `tabix` on the newly-generated VCF file (mazel-tov). We're gonna do something like `ml htslib; tabix -p vcf GTExVCF`. **BUT EVEN BEFORE THAT** I have to gzip our new friend: `bgzip -c GTExVCF.vcf > GTExVCF.vcf.gzip`~~


~~so it looks like this:~~
~~`ml htslib; echo Start Time of Compression; date; bgzip -c GTExVCF.vcf > GTExVCF.vcf.gzip; echo End Time of Compression ; date; tabix -p vcf GTExVCF.vcf.gzip; rm GTExVCF.vcf`~~


~~I could also do `bgzip -f` but I don't want to overwrite the decompressed VCF after *45 minutes* of decryption. *That turned out to be the right call* because I accidentally exited out of the interactive job which had a `screen` that was running the compression. I'm starting over from the dev node. I should have done it on a screen in the dev node but whatever, it's running. *Remind me to run that as a batch script*.~~

~~Okay I cancelled it because it was taking a long time and I don't want it to still be running if I need to leave or anything. Here's what I'm doing: `ml htslib; echo Start Time; date; bgzip -c GTExVCF.vcf > GTExVCF.vcf.gzip; echo End Time; date; tabix -p vcf GTExVCF.vcf.gzip` **Just set it and forget it.** @Wed 12 Dec 2018 04:41:49 PM EST~~

#### Wed 12 Dec 2018 06:57:19 PM EST 
Just found out about `sbatch -d`, will be *extremely* important in writing the master script. S[ecifically, `afterok`. Consult [SLURM documenation](https://slurm.schedmd.com/sbatch.html) for more details.

#### Wed 12 Dec 2018 07:17:42 PM EST 
I just realized that the massive file I was gzipping this whole time was actually the decrypted tar file of what I downloaded. So I got to about 387 GB before I found this out and had to cancel the compression and delete the file. Starting over but instead am doing `tar -xvf phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1.GRU.tar`. Our resulting friend is in a folder already in `.vcf.gz` format, which is what we want.

### 12/11/2018
I found this script from [this stack overflow post](https://stackoverflow.com/questions/2074687/bash-command-to-remove-leading-zeros-from-all-file-names/2074704) that removes all leading zeros from any file in that directory that has a file name that starts with zeros. Not sure exactly how it works and I honestly don't really care. Just hope whoever is using this doesn't use input files that start with zeros.

`shopt -s extglob; for i in *.split; do mv "$i" "${i##*(0)}"; done`

All of the extra code I have to write just to call the LeafCutter scripts is getting out of hand. I really should consider writing a master script soon. Maybe after I call the intron clustering step, I could also call a `rm *.split` to clear the clutter in the working directory, but that's for another time. I should note that the only reason I'm going through all of this trouble is that from all indications, `leafcutter_clustering.py` **only** takes "text file[s] with all junction files to be processed," and not the junction files themselves, unless I'm mistaken which I might be.

Ok, I just checked and it doesn't. Going to `mkdir intronclustering` to reduce working directory clutter. 

~~`python ../clustering/leafcutter_cluster.py -j ${SLURM_ARRAY_TASK_ID}.split -r intronclustering/ -m 50 -o testNE_sQTL -l 500000`~~

From the [documentation](http://davidaknowles.github.io/leafcutter/articles/Usage.html#step-2--intron-clustering)
>This will cluster together the introns fond in the junc files listed in test_juncfiles.txt, requiring 50 split reads supporting each cluster and allowing introns of up to 500kb. The predix testYRIvsEU means the output will be called testYRIvsEU_perind_numers.counts.gz (perind meaning these are the per individual counts).

I also added the `-r` flag to indicate that all of the output files should go to our boy `intronclustering/`. The other flag parameters, I don't know about. For example, I don't know why I should make the minimum cluster reads 50 instead of 30 (which is the default), or why the max intron length should be 50,000bp as opposed to 100,000bp. I'll ask Rajiv if he wants to alter these values but for now I'm going to truck along.

~~Since this is starting to get complicated, I'm going to show what the master script will look like immediately preceeding the above `leafcutter_cluster.py` call:~~


~~`split -d -a 6 -l 1 --additional-suffix '.split' test_juncfiles.txt ''`~~
~~`shopt -s extglob; for i in *.split; do mv "$i" "${i##*(0)}"; done`~~
~~`sbatch intron_cluster.sh`~~


I just realized I'm going to have to do this with GNU-Parallel when I do the full dataset which is frustrating but I guess I gotta if we're going to be efficient. The guy never got back to me about how to best call it but I think I figured it out. Looking at MARCC's documentation, it seems that the most important things are (a) to not use job arrays, (b) `parallel="parallel --delay .2 -j 4 --joblog logs/runtask.log --resume"` and (c) `$parallel "$srun python example_lapack.py {1}000 6" ::: {1..10}`. I feel like meeting with the directors of MARCC made things more complicated for me, but it's okay. I'm actually going to talk to Rajiv about what he thinks about using GNU-Parallel. It might be useful for doing to `.sra` to `.bam` conversion and maybe some of the other more computationally intensive stuff but given that I don't quite understand how it works, I'm not sure how to implement it. **UPDATE:** yeah it's cool I don't need to use it. MARCC is so slow I can't launch a job(s) without first waiting like 4000 minutes.

**Okay, so I just figured out that `leafcutter_cluster.py` doesn't work in a job array, probably because all of the files are merged in the end. Shouldn't be too taxing though, because each `.junc` file is only a few MB in size.**

Ignore what I've striked above. Here is what the code must look like:

```
ml python/2.7-anaconda
mkdir intronclustering/
python ../clustering/leafcutter_cluster.py -j test_juncfiles.txt -r intronclustering/ -m 50 -o testNE_sQTL -l 500000
```

For PCA calculation,
```
python ../../scripts/prepare_phenotype_table.py testNE_sQTL_perind.counts.gz -p 10 # works fine; there's no option for parallelization.
ml htslib; sh testNE_sQTL_perind.counts.gz_prepare.sh
```

Here I am in FastQTL world again. I am going to just download the full genotype matrix in VCF format (~447 Gb): `./prefetch -X 500G --ascp-path '/software/apps/aspera/3.7.2.354/bin/ascp|/software/apps/aspera/3.7.2.354/etc/asperaweb_id_dsa.openssh' ../whole-genome-vcf-cart/cart_prj19186_201812111835.krt > ../../logs/vcf-dl-log.out 2> ../../logs/vcf-dl.err.out` I asked Rajiv which vcf I should download but he's away so I'm just going to do it.

### 12/10/2018
I was able to use `src/12-10-2018/filter_bam.sh` to get rid of all unplaced contigs from the bam files before converting them to `.junc` files and sending them through LeafCutter. `filter_bam.sh` does not remove mitochondrial DNA from the files and I am not sure if LeafCutter can process mitochondiral intron cluters so we will see how this plays out.  

`ls *.filt | tail -n +1 >> filtbamlist.txt`

Doing the same thing; going to round up all of the filtered bam files and put them in this here text file to do the `junc` conversion.

Don't forget to use `wc -l` on `filtbamlist.txt` or `bamlist.txt` or any other text file to be used with a job array to figure out what the job array size should be (if 0 indexed, -1 from size). I need to create some sort of master script that calls all of the job array scripts and basically does everything at once. I made `filter_bam.sh` and `bam2junccall.sh` functional as job arrays today. 

For the intron clustering step, `clustering/leafcutter_python.py -j` intakes a "text file with all junction files to be processed," so I'm trying to figure out how to make this into a job array. Of course, I could make it so each individual line in that text file gets its own text file. 

`split -d -a 6 -l 1 --additional-suffix '.split' <input> ''`

The above command splits a text file. The `-a` flag designates the suffixes (for 5,000 files, 6 should suffice), the `-l` flag indicates lines per output file (1) and the `-d` flag indicates numeric suffixes (`xaa, xab, xac...`) rather than alphabetic suffixes (`x01, x02, x03...`). The `--addtional-suffix` flag ensures that the following command, which would remove leading zeros, can target the resulting files. Lastly, the `''` ensure there is no prefix. The problem with this is that `$SLURM_ARRAY_TASK_ID` generates numbers like `1, 2, 3, etc` but it shouldn't be too difficult to remedy this.

### 12/08/2018
I wrote a script `scratch/filter_bam.txt`, picking up from yesterday. MARCC isn't working right now but I'll try it again tomorrow morning. 

Note: I feel like I didn't mention this before, but `scratch/convertGTEx.txt` is my obviously flawed draft to parallelize converting the GTEx files to `.bam`. I'll have to do the same thing for converting them to `.junc` and pretty much every ensuing step of the project. God what a pain. I'm going to talk to Rajiv about how he thinks we should approach this. 

### 12/07/2018
I met Rajiv about the problem that I wrote about in the README last night without committing the changes and pushing to GitHub (it's okay). Basically the code (`prepare_phenotype.py`) can't handle string inputs that are not simply 'X' or 'Y'. Rajiv slapped together `/data/12-07-2018/GRCh37.bed`, which has the the chromosome number/letter and sequence range (0-terminal) to be used with `samtools` to filter anything that is **not** that from the `bam`. He used the following command: `cat input.txt | sed 's/chr//g' | awk '{print $1"\t0\t"$2}' > GRCh37.bed` and his source was the [UCSC website](http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/referenceSequences/male.hg19.chrom.sizes).

So now I have to figure out which `samtools` commands/flags I need to use to do this. I have to now go back to the previous step, after converting the `sra` files to `bam`s but before turning the `bam` files to `junc`.

From the [samtools documentation](http://www.htslib.org/doc/samtools.html):
>-L FILE
>Only output alignments overlapping the input BED FILE [null].

`samtools view -L GRCh37.bed <file>`

I created a text list of files to parallelize the above command: `ls | tail -n +1 >> bamlist.txt `; I'm going to have `$SLURM_ARRAY_TASK_ID` correspond to each line in the file.

### 12/06/2018
`interact` isn't working on MARCC right now so I'm still using our development node, which should be more than enough in terms of resources. I have converted our 10 `.sras` to `.bams`. Now, following the LeafCutter guide, I'm going to convert them to .junc files, and then I'm going to cluster the introns. 

```
for bamfile in `ls example_geuvadis/*.bam`
do
    echo Converting $bamfile to $bamfile.junc
    sh ../scripts/bam2junc.sh $bamfile $bamfile.junc
    echo $bamfile.junc >> test_juncfiles.txt
done
```

Above is the script provided by LeafCutter to do the .bam -> .junc conversion. Adapted for my needs, the script will look a little something like this:

```
#!/bin/bash
# used in bam2junc.sh
ml samtools
for bamfile in `ls *.bam`
do
    echo Converting $bamfile to $bamfile.junc
    sh /scratch/groups/rmccoy22/leafcutter/scripts/bam2junc.sh $bamfile $bamfile.junc
    echo $bamfile.junc >> test_juncfiles.txt
done
```

bam2junc.sh:

```
#!/bin/bash

leafCutterDir='/scratch/groups/rmccoy22/leafcutter/' ## use only if you don't have scripts folder in your path
bamfile=$1
bedfile=$1.bed
juncfile=$2


if [ ! -z "$leafCutterDir" ]
then
        samtools view $bamfile | python $leafCutterDir/scripts/filter_cs.py | $leafCutterDir/scripts/sam2bed.pl --use-RNA-strand - $bedfile
        $leafCutterDir/scripts/bed2junc.pl $bedfile $juncfile
else
        if ! which filter_cs.py>/dev/null
        then
                echo "ERROR:"
                echo "Add 'scripts' forlder to your path or set leafCutterDir variable in $0"
                exit 1
        fi
        samtools view $bamfile | filter_cs.py | sam2bed.pl --use-RNA-strand - $bedfile
        bed2junc.pl $bedfile $juncfile
fi
rm $bedfile
```

Again, these scripts are altered from what is included in LeafCutter. 

After running this, I got a lot of "problematic spliced reads." Of course, like I said earlier, I'm only using about 10 sra files, so I'm keeping the console output in MARCC. I hope that isn't too big of a deal. The files referenced above are named `callbam2junc.sh` and `bam2junc.sh`, respectively. The first script I pulled from the worked-out example that doesn't actually work that comes with LeafCutter.

I proceeded to the "Intron Clustering" portion of the LeafCutter:

`ml python/2.7-anaconda`

`python /scratch/groups/rmccoy22/leafcutter/clustering/leafcutter_cluster.py -j test_juncfiles.txt -m 50 -o testGTEx -l 500000`

From the [LeafCutter documentation](http://davidaknowles.github.io/leafcutter/articles/Usage.html#step-2--intron-clustering):

>This will cluster together the introns fond in the junc files listed in test_juncfiles.txt, requiring 50 split reads supporting each cluster and allowing introns of up to 500kb. The predix [testGTEx] means the output will be called [testGTEx]_perind_numers.counts.gz (perind meaning these are the per individual counts).

Here's what output `testGTEx_perind_numers.counts.gz` looks like:
```
SRR1068929.sra.bam SRR1068855.sra.bam SRR1069141.sra.bam SRR1068808.sra.bam SRR1068977.sra.bam SRR1069024.sra.bam SRR1069097.sra.bam SRR1069231.sra.bam SRR1069188.sra.bam SRR1068880.sra.bam
14:20813246:20813553:clu_1_NA 20 67 9 30 87 1 2 15 86 4
14:20813285:20813553:clu_1_NA 0 0 0 0 1 0 33 18 0 171
14:20822406:20822968:clu_2_NA 14 48 19 15 55 1 19 29 68 103
14:20822406:20822994:clu_2_NA 5 5 0 1 24 0 6 6 4 7
14:20839791:20840892:clu_3_NA 11 13 32 11 18 1 20 12 1 16
14:20839791:20841170:clu_3_NA 8 3 16 6 30 2 5 10 14 7
14:20841016:20841170:clu_3_NA 7 9 23 7 21 3 19 12 2 8
14:20916970:20917061:clu_4_NA 2 18 5 3 18 0 8 18 115 2
14:20916970:20917123:clu_4_NA 56 93 64 50 102 14 47 53 95 113
14:20917425:20919416:clu_5_NA 2 6 6 1 20 3 1 9 14 0
```

I think each line is an intron and the columns represent the number of split reads from the corresponding file that supports that intron, and the `clu_X` represents which "cluster" it's in. Cool. Now we can get into splicing QTL analysis, which may be the tricky part. I'm doing this all with just 10 SRA files and it's really highlighting the need to parallelize this process, because if we don't, it would take like 20 years just to get to this part using the full dataset.

[Splicing QTL Analysis](http://davidaknowles.github.io/leafcutter/articles/sQTL.html):

`python /scratch/groups/rmccoy22/leafcutter/scripts/prepare_phenotype_table.py /scratch/groups/rmccoy22/bamfiles/testGTEx_perind.counts.gz -p 10`

I get this error message:
```
[aseyedi2@jhu.edu@rmccoy22-dev bamfiles]$ python /scratch/groups/rmccoy22/leafcutter/scripts/prepare_phenotype_table.py testGTEx_perind.counts.gz -p 10
Starting...
Parsed 1000 introns...
Traceback (most recent call last):
  File "/scratch/groups/rmccoy22/leafcutter/scripts/prepare_phenotype_table.py", line 171, in <module>
    main(args[0], int(options.npcs) )
  File "/scratch/groups/rmccoy22/leafcutter/scripts/prepare_phenotype_table.py", line 99, in main
    chrom_int = int(chr_)
ValueError: invalid literal for int() with base 10: 'GL000219.1'
```
Going to tell Rajiv about it.

### 12/05/2018
Got up this morning and looks like firing off the command `for f in $PWD/*.sra; do ./sam-dump $f | samtools view -bS - > $f.bam; done` with a 4-hour long interactive session wasn't even enough for just 10 files. I'm converting the rest of the files on `rmccoy22-dev`. 

### 12/04/2018
We met with the directors of MARCC last week and we talked about using a tool called `GNU-Parallel` which, from what I remember, is a good way to multithread jobs in a way such that if one crashes, we can pick up from where it left off or at least figure out what the problem was. Also, Keven wanted us to measure for him the download speed of the SRA files so he has an idea of how fast a TB downloads. I'm gonna do a dry run with LeafCutter and only about 10 files to see how it works, first. Then I'm going to study `GNU-Parallel` a bit and then email Kevin saying I don't understand it. For the purposes of this practice run, I'm using `for f in $PWD/*.sra; do ./sam-dump $f | samtools view -bS - > $f.bam; done` to convert the boys into `.bam` files.

## Update
`for f in $PWD/*.sra; do ./sam-dump $f | samtools view -bS - | tee $f.bam; done`

I tried out the above command and it was a nightmare. I thought it would be readable for whatever reason. It would be cool if I figure out how to add like a progress bar or something just so I can monitor the progress of the conversion, but I guess I won't even have to with parallelization.

### 11/27/2018
Using `nohup` didn't work, or something else happened, but the cart files did not download. Disappointingly. I'm going to try again, but this time with `screen`. I'm going to use `./prefetch -X 50000000000 /scratch/groups/rmccoy22/sratools/cart_prj19186_201811221749.krt > /scratch/groups/rmccoy22/logs/SRA_DL.out`, with the output redirected to a `logs` folder so that the download progress can be monitored by `tail -f`. Just to clarify, I am using the dtn node. Also, Ubuntu's built-in text editor isn't working on my laptop, weirdly.

### 11/26/2018
The more things change, the more they stay the same. We are meeting with the director of MARCC on Wednesday to discuss our options when dealing with big ol' datasets. Looks like using `~/scratch` is a no-go, but `~/work` apparently has up to 50TB of space AND is a scratch partition, so that's what we're going to use. We're still waiting to meet with the dude on Wednesday but in the meantime, I'm going to play around with a small subset of the total data to see what's what. I already made some progress last week; I came up with a little ditty like this: `for f in $PWD/*.sra; do ./sam-dump $f | samtools view -bS > $f.bam; done`. However, this one-off command doesn't utilize parallelization, so Rajiv recommended something that makes use of SLURM's capabilities:
```
#!/bin/bash
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

echo "Finished with job $SLURM_ARRAY_TASK_ID"
```
Since we won't need to delete and move things back and forth, I have commented-out the parts of the code that are now both practically and conceptually obsolete. I'm going to use the following command: `nohup ./prefetch -X 50000000000 /scratch/groups/rmccoy22/sratools/cart_prj19186_201811221749.krt &` to download the boys (the boys being the SRA files).

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

