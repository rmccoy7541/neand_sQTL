## Updates

These updates are read from most recent date at the top to initial entry at the bottom.

REMEMBER
- You need to make the pipeline generalizable, especially the part after separating the phenotype table by tissue.

### 9:23 AM 4/5/2019
Holy crap... so splitting the TESTIS had a mistake, in that many of the files were split along the middle of the line, destroying the whole thing. I'm going to have to split again and then do NPEC, but this time split by lines and not by bytes.

I fucked up and accidentally deleted the TESTIS nominal pass concatted file. Damn. Whatever, it was becoming a pain in the ass to do it the hard way. I'm just going do do a NomPass and then run NPEC and get the Neaderthal-extracted sQTLs the normal way, which is from the separate outputs. Whatever.

Bamming the WHLBLD sra's

### 4:29 PM 4/4/2019
I'm trying to call NPEC on all of the remaining tissues (BRNCHA, TESTIS, WHLBLD). I'm going to have to do all of WHLBLD all over again because I don't have the intermediary files. Which is fine for now since my focus is putting together a presentation. For TESTIS, the nominal was close to 1 terabyte, so I broke it up into ~ one hundred pieces and am going to try to run NPEC on that. I think BRNCHA is (almost) done with NPEC so I'm going to cat them and put them in `~/data`. So really I have one tissue to work on immediately (TESTIS). I think the plan is talk about some of the top Neanderthal genes by tissue. 

### 04/02/2019
I finished with the three tissues using the SPrime calls (BRNCTXB, LIVER, and MSCLSK), and MSCLSK has 21 sQTL hits, whereas the others have 0. 

```
for i in {1..100}; do
   cat BRNCTXB_nominals_chunk_${i}_out.txt | gzip -c >> BRNCTXB.nominals.all.chunks.NE_only.txt.gz
done

sbatch --export=listPath=$PWD,tissue=$(echo BRNCTXB),scripts=$scripts ../NomPassExtractCall.sh

Rscript ${scripts}/R/QQPlot-Viz.R $PWD res_SPRIME/BRNCTXB.nominals.all.chunks.NE_only.txt.gz BRNCTXB.permutations_full.txt.gz ${data}/../analysis/SPRIME/sprime_calls.txt
```

### 04/01/2019
Happy April fools day. I (speaking on behalf of MARCC) finished producing the nominals again and I'm goint to try concatting (again again) to make sure it was all good. The MSCLSK and BRNCTXB finished concatting, I'm going to `sha1sum` them to see if they match the ones I already made. I don't know.

Okay, surprise, none of the other tissues finished concatting. What the hell. I feel like I'm wasting my time.
```
for i in {1..100}; do
   cat LIVER_nominals_chunk_${i}.txt | gzip -c >> LIVER.nominals.all.chunks.in-screen.txt.gz
done
```
```
for i in {1..100}; do
   cat BRNCTXB_nominals_chunk_${i}.txt | gzip -c >> BRNCTXB.nominals.all.chunks.in-screen.txt.gz
done
```
I'm concatting again, but in a supervised way this time. Hopefully 6 hours and 4 cores is enough to do this. Not having the dev node is really screwing things up.

I ran NPEC on all of the generated nominal chunks for each one of the three tissues I'm working with. I may eventually need to do it for the other tissues I've already worked with but we will see. The goal is to produce NPE files to then concatenate and extract the top Neanderthal calls from. I'm doing practice-practice presentation next Monday so I really need to step it up. I don't know why I volunteered for it.


### 03/31/2019
```
for i in {1..100}; do
   cat MSCLSK_nominals_chunk_${i}.txt | gzip -c >> MSCLSK.nominals.all.chunks.in-screen.txt.gz
done
```
Okay something fucked up. I'm not even sure what I'm doing anymore. I got this weird message about some of the files in the LIVER directory not being accessible so I impulsively deleted them. Not sure that I had to do that but I did it. 

### 03/30/19
I just finished generating the nom pass for LIVER and concatting them and the new concatted file is different from the old one. However, I concatted them on a screen that crashed, so I'm not sure if it went all the way through. I'm going to concat them again and see if they're any different. Important to note that the old and new LIVER nominals concatted file are about the same size but have a different shasum. 

### 03/29/19
I feel like there's something fishy about the two NE only nominal passes (BRNCTXB and LIVER) I just generated... They're both less than 2 Mb, whereas the other tissues' NE only nom pass are like... at least 100 Mb. There is something fishy about this. I feel like I should redo QTLtools steps for BRNCTXB, LIVER, and MSCLSK manually. 

I'm going to try generating NE only files for MSCLSK. If that doesn't work, I'm going to redo the QTLtools step for all three of them interactively (I've already requested a 5-day-long interactive session) and worry about debugging the master script later. 

Update: the MSCLSK output is only 12 Mb. That's not good. Once my unlimited interactive job gets approved, I'm going to try it again, and by 'it' I mean running QTLtools but interactively this time.

### 03/28/19
I just talked to Ryan Bradley over at his office at the CS building. He basically told me that I should use `snakemake`. Also that I could get away with calling an unlimited `sbatch` and then maybe using sbatch to do the rest of it. Make it just a series of `sbatch` calls. That could work. 

The dev node is down and I still need to do some stuff with the 6 tissues we have right now. I'm "focusing" on the presentation for RECOMB-GEN so I don't really have all that much time to spend on optimizing the pipeline.

### 12:23 PM 3/22/2019
Dealt with concatting the remaining permutation pass files that didn't work for MSCLSK. Calling another conditional pass. Also doing nom pass extract for MSCLSK, BRNCTXB, and LIVER. All I have left with these files is QQviz and I'm done. Though it won't matter since I'm just going to do this all over again with the fully functioning pipeline.

Speaking of which, the lung, thyroid and skin files are still downloading. Once they're done, I'm going to call the pipeline on them and see what happens.


### Thu 21 Mar 2019 03:21:19 PM EDT 
Fixed conditional pass and nompassextract and now they're about to run. My next step is to download all of the files for lung, thyroid and exposed skin and try to pipeline on them since the quota is about to run out. 

```
{
for i in $(ls *krt)
do 
   /scratch/groups/rmccoy22/progs/sra-tools/bin/prefetch -O "/home-1/aseyedi2@jhu.edu/work/Ne_sQTL/sra/lung_skinEx_thy" -X 500G --ascp-path '/software/apps/aspera/3.7.2.354/bin/ascp|/software/apps/aspera/3.7.2.354/etc/asperaweb_id_dsa.openssh' $i
done
} > DL-prog.out 2> DL-prog.err
```
that should work. Going to change the sra variable to this location.

### Wed 20 Mar 2019 03:39:04 PM EDT 
Conditional pass and NomPassExtract are failing because of some errors in the code that I have since fixed. 

### Tue 19 Mar 2019 09:14:28 AM EDT 
Running the QTLtools for-loop on the `unlimited` partition, let's see how this goes.

### Mon 18 Mar 2019 12:23:25 PM EDT 
Everything is coming along very well. I'm making the pipeline generalizeable with minimum intervention. Here are the lines of code I used to try to make the covariate files useable:

`for i in GTEx_Analysis_v7_eQTL_covariates/*; do echo $i | awk -F'[./]' '{print $2}' >> GTExCovNames.txt; done`

For some reason, this for loop doesn't work:
```
for line in $(cat tissuesused.txt)
do
   head -1 $line/1_$line.txt > $line/$line.phen_fastqtl.bed
   echo "Concatenating $line phenotypes..."
   for file in $(ls $line/*)
   do
      cat $file | sed -e1,1d >> $line/$line.phen_fastqtl.bed
   done
done
```

What happens is that I end up with an enormous, 700Gb+ file in BRNCTXB. Each one of the files is a few dozen Mb. I feel like there's a recursion problem happening but I can't figure it out.
```
for line in $(cat tissuesused.txt)
do
   head -1 $line/1_$line.txt > $line/$line.phen_fastqtl.bed
   echo "$line/1_$line.txt > $line/$line.phen_fastqtl.bed"
   for file in $(ls $line/*)
   do
      echo "$file goes into $line/$line.phen_fastqtl.bed"
      cat $file | sed -e1,1d | head
   done
done
```

I figured it out. The for-loop would recur on the concatted file because of `ls $line/*`. Should be `ls $line/*_*.txt`.

Okay, I was able to get up to the QTLtools part of things. I'm almost finished.

I'm running into a problem now:

I can run the QTLtool part with little problem except for this: blah blah blah something about a token near `done`.


### Sun 17 Mar 2019 06:45:59 PM EDT 
Gonna do bam -> junc. I did it once and it worked but I forgot to up the number from 318 to 734, so I'm going to do it again.

### 03/16/2019
#### Sat 16 Mar 2019 06:16:56 PM EDT 
Just finished bamming the muscle files. `samtools quickcheck` reveals three bams that are corrupt:
```
[aseyedi2@jhu.edu@compute0036 skeletal_muscle]$ samtools quickcheck *bam
SRR2135308.sra.bam had no targets in header.
SRR2135331.sra.bam had no targets in header.
SRR2135378.sra.bam had no targets in header.
```
There's probably just something wrong with those files. The sizes for those sra's are only a few dozen MB. I removed the bams and put all of the bams in a single folder to continue with (as I intended).

After moving them all into one folder, I did another `quickcheck` just to be sure. Some of the files that I tried to reprocess are still there and broken. I'm going to remove them.

```
[aseyedi2@jhu.edu@compute0036 frontallobe_liver_muscle]$ cat quickcheckfail.txt | cut -d' ' -f1 | cut -d'.' -f1,2,3 > failedbamconversions.txt 
[aseyedi2@jhu.edu@compute0036 frontallobe_liver_muscle]$ cat failedbamconversions.txt 
SRR2135317.sra.bam
SRR2135333.sra.bam
SRR2135349.sra.bam
SRR2135358.sra.bam
SRR2135365.sra.bam
SRR2135366.sra.bam
SRR2135380.sra.bam
SRR2135387.sra.bam

```

Between frontal lobe, skeletal muscle and liver tissues, there are 734 samples. 

### 03/15/2019
#### Fri 15 Mar 2019 11:25:03 AM EDT 
The sra's were successfully converted to bam files for the liver and frontal cortex samples. Skeletal muscle samples finished downloading, 479 of them.

The following sra's failed in conversion:
```
SRR1384049.sra.bam was missing EOF block when one should be present.
SRR2135317.sra.bam had no targets in header.
SRR2135333.sra.bam had no targets in header.
SRR2135349.sra.bam had no targets in header.
SRR2135358.sra.bam had no targets in header.
SRR2135365.sra.bam had no targets in header.
SRR2135366.sra.bam had no targets in header.
SRR2135380.sra.bam had no targets in header.
SRR2135387.sra.bam had no targets in header.
```

`samtools quickcheck *bam 2> quickcheckfail.txt`
`cat quickcheckfail.txt | cut -d' ' -f1 | cut -d'.' -f1,2 > failedbamconversions.txt`

Still getting the same error. These files are just broken or something. Doing the bam conversions for skeletal muscle right now. 

### 03/13/2019
#### Wed 13 Mar 2019 12:59:37 PM EDT 
I'm going to start with liver and frontal lobe tissues. I downloaded a cart file with 144 liver files and 122 frontal cortex samples. However, that's less than the supposed amount according to sample attributes DS. I'm going to talk to Rajiv about this.

### 03/12/2019
#### Tue 12 Mar 2019 12:03:20 PM EDT 
My computer crashed but basically the BRNCHA and TESTIS are done and am gathering analysis now.

### 03/11/2019
#### Mon 11 Mar 2019 10:58:20 AM EDT 
I forgot to write about it but basically I finished the permutation step for both tissues. Finished concatting for BRNCHA, now gotta do the same for TESTIS. Then I'm going to pull the FDR for BRNCHA.

```
sbatch --export=VCF=$VCF,pheno=$pheno,tissue=$(echo TESTIS),covariates=$(echo Testis.v7.covariates_output.txt),permutations=$(echo TESTIS.permutations_full_FDR.thresholds.txt) ${scripts}/sh/CondPass.sh

sbatch --export=VCF=$VCF,pheno=$(echo BRNCHA.pheno.bed.gz),tissue=$(echo BRNCHA),covariates=$(echo Brain_Cerebellum.v7.covariates_output.txt),permutations=$(echo BRNCHA.permutations_full_FDR.thresholds.txt) ${scripts}/sh/CondPass.sh
```


### 03/07/2019
#### Thu 07 Mar 2019 02:02:05 PM EST 
I just finished concatting and compressing the nominals for TESTIS and BRNCHA. They're both at around ~120Gb, which is not bad but still big for the number of samples. I'm going to perform a sanity check aka I'm going to try to reproduce nominal chunk 50 for both tissues and if the `sha1sum` doesn't match up, I'm going to retrace my steps.

`/scratch/groups/rmccoy22/progs/QTLtools/QTLtools_1.1_Ubuntu14.04_x86_64 cis --vcf $VCF --bed "../${pheno}" --cov  "../Brain_Cerebellum.v7.covariates_output.txt" --nominal 1 --chunk 50 100 --out "BRNCHA_nominals_chunk_50_SANITYCHECK.txt"`

BRNCHA is good

```
[aseyedi2@jhu.edu@rmccoy22-dev sanitycheck]$ sha1sum BRNCHA_nominals_chunk_50_SANITYCHECK.txt 
c47b1b6868b59868bfee98c29736aa618439b382  BRNCHA_nominals_chunk_50_SANITYCHECK.txt
[aseyedi2@jhu.edu@rmccoy22-dev sanitycheck]$ sha1sum ../BRNCHA_nominals_chunk_50.txt 
c47b1b6868b59868bfee98c29736aa618439b382  ../BRNCHA_nominals_chunk_50.txt
[aseyedi2@jhu.edu@rmccoy22-dev sanitycheck]$ cd ../../TESTIS/

```
Going to try it again for TESTIS
`/scratch/groups/rmccoy22/progs/QTLtools/QTLtools_1.1_Ubuntu14.04_x86_64 cis --vcf $VCF --bed "${pheno}" --cov  "../Testis.v7.covariates_output.txt" --nominal 1 --chunk 50 100 --out "TESTIS_nominals_chunk_50_SANITYCHECK.txt"`
```
[aseyedi2@jhu.edu@compute0003 sanitycheck]$ sha1sum TESTIS_nominals_chunk_50_SANITYCHECK.txt 
821da43b75f23fe802f70d8521e618bf94ffd4f1  TESTIS_nominals_chunk_50_SANITYCHECK.txt
[aseyedi2@jhu.edu@compute0003 sanitycheck]$ sha1sum ../TESTIS_nominals_chunk_50.txt 
821da43b75f23fe802f70d8521e618bf94ffd4f1  ../TESTIS_nominals_chunk_50.txt

```

Okay, TESTIS is good too. There is no obvious problem. I'm going to call the nominal pass extract for BRNCHA now.

Both of the NomPassExtractCalls are running. Now we wait.

### 03/06/2019
#### Wed 06 Mar 2019 06:07:10 PM EST 
I'm concatenating the QTLtools nominal outputs for testis and cerebellum and they're coming out to be massive files compared to whole blood. I'm going to have to look into this. For contrast, the total size for the concatenating nominals for WHLBLD were 54G. 

Testis just stopped at around 130Gb, BRNCHA is still going at ~80Gb. I'm running Testis again just to be sure.

Double the size isn't that bad though.

I still absolutely need to find out how many files were included for each tissue in the analysis freeze.


### 03/05/2019
#### Tue 05 Mar 2019 11:15:07 AM EST 
Working on `sraNameChangeSort.sh`.

Done.

I need to figure out how to make the part all after separating the samples by tissue generalizeable. Right now it's all just being done manually in the interest of time. But this is something I need to remember to do.

QTLtools nominal pass is running right now on both the brain - cerebellum as well as testis files on the dev node. It's very important that I remember to figure out how to make this generalized because right now I'm doing it manually. Going to check back in a few hours and do the permutation pass, which also needs to be generalized. Maybe some sort of user input?

### 03/04/2019
#### Mon 04 Mar 2019 06:38:28 PM EST
Forgot to write in here all day but it's okay. `sraNameChangeSort.sh` doesn't work; I feed it only cerebellum and testes tissue samples and it gives me samples matching to every tissue, which isn't okay.

### 03/01/2019
#### Fri 01 Mar 2019 10:56:52 AM EST 
I'm going to try to qqplot the results. Rajiv gave me a script I could use in service of that mission.

`cat GTEx_v7_Annotations_SampleAttributesDS.txt | awk -F '\t' '{if ($17 == "RNASEQ" && ($7 == "Brain - Cerebellum" || $7 == "Testis")) print $1, $7}' | sort > gtex_cerebellum_testis.txt`

`join -1 2 -2 1 dbgap_rnaseq_manifest.txt gtex_cerebellum_testis.txt | awk '{$7=$8=$9=""; print $0}' | awk '{gsub("Testis", "", $5); print $0}'`

Rajiv made the file that contains all of the rna-seq samples from the dbgap manifest. I used it to extract SRRs and am right now converting an SRA that was not converting initially. We're just plowing ahead, and since the corrupted files were not in the analysis freeze, we wouldn't have included them anyway.

SRR1358260.sra.bam doesn't work

### 02/28/2019
#### Thu 28 Feb 2019 11:36:08 AM EST 
I'm not sure if I can really do anything when it comes to the cerebellum and testis files. Again, I'm sure I could just find out which among those files were and were not included in the v7 analysis freeze, and if none of them are the corrupted files, I could move ahead with the process.

I'm doing conditional pass for the whole blood samples now.

Done. I have the results for the whole blood samples.

### 02/27/2019
#### Wed 27 Feb 2019 10:41:07 AM EST 
All of the SLURM jobs are done. 

For the permutation pass, I got an error - minor stuff. Fixed it and am trying it again.

I need to figure out how to streamline the concatenating of the covariates part.

As for the testis files, I got several errors and 5 broken files according to `samtools quickcheck`:
```
[aseyedi2@jhu.edu@compute0195 testis]$ samtools quickcheck *bam
SRR2135285.sra.bam had no targets in header.
SRR2135301.sra.bam had no targets in header.
SRR2135302.sra.bam had no targets in header.
SRR2135320.sra.bam had no targets in header.
SRR2135354.sra.bam had no targets in header.
SRR2135377.sra.bam had no targets in header.
```

I just tried converting these files independently and still same error. Going to delete them and include a quickcheck step in the master script.

Rajiv just told me to try my best to have the same set of files used in the analysis v7 freeze, so I will try my dangest to get those files.

Okay, I've tried redownloading those files and the problem is that they're each only a few dozen MB in size and I get the same file each time I download. So this is a problem. I'm going to email NCBI.

Meanwhile, I got the analysis freeze for the cerebellum and testis files: 

`wget https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SampleAttributesDS.txt`
`cat GTEx_v7_Annotations_SampleAttributesDS.txt | awk -F '\t' '{if ($17 == "RNASEQ" && $7 == "Brain - Cerebellum" || $7 == "Testis") print $1, $7}' | sort > gtex_cerebellum_testis.txt`

I probably have to include a seperate part for extracting the tissues of interest.

I'm going to try redownloading the files for both cerebellum and testis and see if that helps with the problem. Didn't work.

I emailed NLM about this problem so I guess I'll have to wait until they get back to me. I guess I could continue with the pipeline and hope the files are not in the analysis freeze downstream. PermPass is still running.

### 02/26/2019
#### Tue 26 Feb 2019 11:23:16 AM EST 
I'm still having major difficulty doing the sra2bam conversion in a stream-lined way, so I'm just trying that over and over again hoping something will eventually work. 

The WHLBLD nom pass extract did NOT work, no package called 'tidyverse' also some other errors that I don't know what to do about. I tried fixing it by loading `gcc`. Let's see if that works.

Meanwhile, for the testis files, the conversion is finally working, but `SRR1069734.sra` is having trouble converting. I'm getting the following error:

```
...
2019-02-26T17:10:12 sam-dump.2.9.2 err: unknown while creating file within file system module - unknown system error 'Cannot send after transport endpoint shutdown (108)'
2019-02-26T17:10:12 sam-dump.2.9.2 err: unknown while creating file within file system module - error with https open 'https://sra-download.ncbi.nlm.nih.gov/traces/refseq/NC_000011.9'
2019-02-26T17:10:12 sam-dump.2.9.2 int: unknown while creating file within file system module - VCursorCellDataDirect( row#1921079 . idx#3 . READ ) char_ptr failed
[W::sam_read1] Parse error at line 39693787
[main_samview] truncated file.
```

Not sure what that means because I don't recall there being any problem with downloading any of the testis files. I'm going to try downloading that file by itself and then running a conversion on the dev node.

I think `NomPassExtract` worked... But got some OOM errors. None of the memory flags work for SLURM anymore so I'm going to resubmit the job with multiple cores per array task. Fingers crossed. Okay, no matter what I do, I can't seem to avoid an OOM error for `NomPassExtract`. I sent an email to the MARCC help desk though I don't anticipate them being much help. Let's see what happens.

Rajiv downloaded that one file `SRR1069734.sra` for me and now I'm trying to do the sra2bam conversion independently on the dev node. Let's see how this will pan out. I'm sure there are other sra's that didn't fully convert.

The `NomPassExtract` step is almost done; one file did not finish in the 1 hour of time that was allotted. I'm running chunk 32 separately [sp?] on the dev node screen.

Strange; it took over an hour on SLURM but was over super quick on the dev node.

### 02/25/2019
#### Mon 25 Feb 2019 11:32:00 AM EST 
WHLBLD is done with its nominal pass. Concatting now.

The testis bam conversion still didn't work.

I figured out what the problem was: I misnamed the directories to be exported into the script, which is a shame because this script takes forever.

Still concatting WHLBLD nominals. Just finished. Now I'm trying to extract the Neanderthal sequences. 

Okay, we're trying the testis sra -> bam conversion once again as well as the nominal pass extract call. We'll see what happens; MARCC has been frustratingly slow today.

#### Mon 25 Feb 2019 05:51:55 PM EST 
`NomPassExtract.R` didn't work because of some technicality. I fixed it, I hope. I'm just trying to generalize this pipeline as much as possible so that it requires the least amount of manual input but it's not quite that easy.


### 02/21/2019
#### Thu 21 Feb 2019 12:56:35 PM EST 
I haven't been great at keeping up with the log lately, but basically we're redoing the WHLBLD pipeline but filtering out NO introns for `prepare_phenotype_table.py`. 

I just submitted `QTLtools-Filter.sh` on the WHLBLD samples. While I have that going on, in parallel, I want to convert - Okay, the sra -> bam conversions for the cerebellum sra's were all broken pretty much:

```
[aseyedi2@jhu.edu@bc-login01 brain_cerebellum]$ egrep -L -I "err|fail" slurm-32611501_* | wc -l
129
```

Moving testis over to the project directory. Shifting gears back to WHLBLD.

Two files did not survive the conversion to bam.
```
SRR2135355.sra.bam had no targets in header.
SRR2135374.sra.bam had no targets in header.

```
I'm just going to get rid of them since I've already tried to make them work.

Running QTLtools nominal pass on the full set of introns for whole blood. 

Okay, I'm going to head out. Cerebellum files have been converted except for those two files, testis files are now converting to bam, and we're doing a nominal pass for the whole blood but with all introns ever. I'm going to throw the cerebellum and testis files together after this step.

### 02/19/2019
#### Tue 19 Feb 2019 02:55:46 PM EST 
We did it. We went from beginning to end on one set of tissues (whole blood). I'm going to do this again but for brain tissues and testis, and I'm going to split them up by brain region.

UPDATE: Okay dbGaP is broken. For the time being, I'm going to rerun leafcutter but with a threshold of 0% for the intron excision filter. Okay, I'm running `sraNameChangeSort.R` right now but on a ton of different introns.

### 02/14/2019
#### Thu 14 Feb 2019 08:14:27 AM EST 
I forced normal distro for perm pass and that apparently fixed the problem. I'm now re-running conditional pass.

Conditional pass done. Finally, the results for this project have been generated. Time for the analysis.

#### Thu 14 Feb 2019 01:58:55 PM EST 
I want to experiment with keeping all of the scripts I use in the pipeline just floating around the `/src` directory and see if that makes writing the pipeline any easier.

### 02/13/2019
#### Wed 13 Feb 2019 09:36:34 AM EST 
Produced a permutation pass. Doing FDR correction:

`Rscript ../../../../../progs/QTLtools/script/runFDR_cis.R permutations_full.txt.gz 0.05 permuatations_full_FDR`

```
Processing fastQTL output
  * Input  = [ permutations_full.txt.gz ]
  * FDR    =  0.05 
  * Output = [ permuatations_full_FDR ]

Read Input data
  * Number of molecular phenotypes = 99351 
  * Number of NA lines = 0 
  * Correlation between Beta approx. and Empirical p-values = 0.9703 

Process Input data with Qvalue
  * Proportion of significant phenotypes = 0 %

Determine significance thresholds
  * Corrected p-value threshold =  0.00494386 
There were 50 or more warnings (use warnings() to see the first 50)
There were 50 or more warnings (use warnings() to see the first 50)
  * pval0 =  NaN  +/-  NA 
  * test0 =  NaN  +/-  NA 
  * corr0 =  NaN  +/-  NA 
  * test1 =  NaN  +/-  NA 
  * pval1 =  NaN  +/-  NA 

Write significant hits in [ permuatations_full_FDR.significant.txt ]

Write nominal thresholds in [ permuatations_full_FDR.thresholds.txt ]
```

#### Wed 13 Feb 2019 10:37:04 AM EST 
According to Rajiv, we have a problematic distribution of p-values. However, LC forces QQ normalization for its phenotype outputs, so I have to force normal distribution for the permutation and conditional passes.


### 02/12/2019
#### Tue 12 Feb 2019 09:39:00 AM EST 
The perm pass is taking way too dang long. I'm running it as a batch script now. Rajiv brought to my attention that `NomPassExtract.R` only produced 97 output files, whereas I need 100.

biomart links to ensembl annotation to pull gene info

what are the #s LC puts in phen table

`PermPass.sh` did not have enough time. I upped it from 3 hrs to 12 hrs (just in case). Below is how I was able to extract the numbers for the chunks whose permutation pass failed.

`sacct | grep "TIMEOUT" | awk '{ print $1 }' | cut -f2 -d _ > failedpermpass.txt`

`awk -vORS=, '{ print $1 }' failedpermpass.txt | sed 's/,$/\n/'`

`1,2,4,8,10,13,14,19,20,21,22,23,28,29,30,31,32,33,34,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,69,83,88`

#### Tue 12 Feb 2019 06:33:18 PM EST 
I messed up `PermPass.sh`. Fixed it, check out version history.


### 02/11/2019
#### Mon 11 Feb 2019 03:14:59 PM EST 
Did not document on Friday so not sure what happened. We are now concatenating the nominal pass files while Rajiv performs some manipulations in R to try to figure out how to pull Neandertal sQTLs from the final product.

UCSC Genome Browser hg19

GGV - geography of genetic variants

mpi archaic genome browser

#### Mon 11 Feb 2019 06:35:18 PM EST 
We covered a lot of ground today. I'm going to go ahead and process the QTL nominal pass results with the script he gave me and go ahead with the proceeding steps of QTLtools.


#### Mon 11 Feb 2019 08:49:17 PM EST 
Made some progress. Started on the permutation pass, filtered out the Neanderthal sequences from the nominal pass. Going to take a break now.


### 02/07/2019
#### Thu 07 Feb 2019 10:47:43 AM EST 
Rajiv got the full list of SRAs that were used in the analysis including the ones with ambiguous IDs, so I have to figure out which ones they are (he sent me a table via slack) and then run LC again. Almost done. Don't forget. I also have to figure out how to change the headers on the final phenotype files for LC to the subject ID, but that's down the line.

Apparently, none of the SRRs he sent me via slack in `whole_blood_analysis_freeze_genotyped.txt` match up with the ones I already have, which is a bit fishy. I messaged Rajiv about it, we'll see what he says.

Okay the above road is a dead-end. We made a new file, `genotyped_samples_used.txt`. I'm using the SRRs from that because those include genotyped files. We are back to the step immediately following junc generation.

Using `genotyped_samples.txt`, n = 355

Okay, I'm at the point after getting `tissue_table.txt` but before `sraNameChangeSort.R`. 

I found a file named just `.txt`. I don't know what it's from, but I have a hunch that it's from `sraNameChangeSort.R`. It's not blank. I will have to investigate this.

Nevermind, I forgot that I'm changing the COVARIATE file.

Okay I submitted `QTLtools-NomPass.sh`. I might finally be on the forefront of this project again.

`zcat nominals.chunk*.txt.gz | gzip -c > nominals.all.chunks.txt.gz` to concatenate the outputs.

### 02/06/2019
#### Wed 06 Feb 2019 07:37:21 AM EST 
Our new friends are done converting. I'm going to catch them up along the pipeline.

Filtering now.

SRR5125340
SRR5125397
SRR5125400
SRR5125580
SRR5125595

I'm right before the phenotype mapping step.

SRR607214 got left behind; need to convert to bam, filt, junc.

`./sam-dump SRR607214.sra | ./samtools view --threads 23 -bS - > SRR607214.sra.bam`

Once I catch this last remaining file up to speed, I am ready to prepare phenotype table.

#### Wed 06 Feb 2019 11:02:00 AM EST 
Still processing.

#### Wed 06 Feb 2019 12:13:38 PM EST 
Done, and up to intron clustering step. Now I have to use `sorted_SRRsNeeded.txt` (n = 391) to do junc stuff.

`python ../../../aseyedi2/leafcutter/clustering/leafcutter_cluster.py -j sorted_SRRsNeeded.txt -r intronclustering/ -m 50 -o Ne-sQTL -l 500000`

Used this line to concatenate the whole blood samples:

`for q in {1..22}; do echo "Chr $q..."; awk 'FNR==1 && NR!=1{next;}{print}' "${q}_WHLBLD.txt" >> WHLBLD.txt; done;`

I did a bunch of stuff that I didn't log because it was a pain the the ass, but basically nothing out of the ordinary, I just followed the steps in the master script but applied them to just one tissue type. Now I'm running into a problem where I run `mergePCs.R` with both the GTEx-supplied covariates and the LC-generated covariates as inputs and I'm getting an essentially blank file with all of the appropriate headers as the output. I sent the inputs and the script over to Rajiv to get his insight. 

#### 3:22 PM 2/6/2019
Okay the problem is that I use different ID schemes in the two PC files, so the GTEx one is in subject ID whereas the LC one is in SRR. pretty obvious in hindsight. I have to convert the headers to the subject ID to match the GTEx scheme.

-------

Below are all the commands we ran yesterday. We did too much for me to keep perfect track of. I only included the most important ones.
```
cat GTEx_v7_Annotations_SampleAttributesDS.txt | awk -F '\t' '{if ($17 == "RNASEQ" && $7 == "Whole Blood") print $1}' | sort > whole_blood_analysis_freeze.txt

cat WholeBlood.txt | cut -f21,23 | sed -e 1,1d | sort -k2,2 > srr2gtex_sampleid.txt

join -1 2 -2 1 srr2gtex_sampleid.txt whole_blood_analysis_freeze.txt | awk '{print $1"\t"$2}' > matchedIDs.txt

matchedIDs.txt | tr '-' '\t' | awk '{print $1"-"$2"\t"$1"-"$2"-"$3"-"$4"-"$5"\t"$6}' > subject_sample_srr.txt

zcat GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz | head -1000 | grep '#' | tail -1 | cut -f10- | tr '\t' '\n' | sort > genotyped_subjects.txt

cat subject_sample_srr.txt | awk -F"\t" '!_[$1]++' > samples_used.txt
```

### 02/05/2019
#### Tue 05 Feb 2019 10:55:40 AM EST 
I was working on my windows laptop last night and forgot to push the changes to GitHub, but I basically remember what I was stuck on. I have a file named `SRRs.sorted` which contains all of the SRR IDs with the GTEx sample IDs for all whole blood samples. What I want to do is compare those to the sample IDs that were used for eQTL analysis, get the matching GTEx IDs and turn those back into SRR IDs to then determine which of those I've already converted and if I need to convert any other files to `.junc`.

Okay, so the analysis freeze is 11,688 samples big. I now have to extract all of the samples that match the sample IDs of the whole bloods I'm working with.

`awk 'FNR==NR {a[$1]=$2; print a[$1]; next}; $1 in a {print a[$1]}' SRRGTEX.txt SRRs.txt`
`919`

Why 919? That's too many. There should be fewer than that.

#### Tue 05 Feb 2019 12:40:51 PM EST 
Rajiv and I have determined that there are some files that we cannot access that were used in the eQTL analysis. There are also several SRAs that correspond to some of the same samples.

```
GTEX-V955-0005-SM-3P5ZC	SRR5125400
GTEX-V955-0005-SM-3P5ZC	SRR5125401
GTEX-V955-0005-SM-3P5ZC	SRR810595
```

These guys are all duplicates.

There is a new file called `samples_used.txt` in the `data/` directory. There, you will find our map for samples IDs/SRRs/subject IDs.

```
IFS=$'\n'       # make newlines the only separator
for sample in $(cat samples_used.txt)
do
    echo $sample | cut -f3
done
```

#### Tue 05 Feb 2019 01:29:15 PM EST 

```
[aseyedi2@jhu.edu@rmccoy22-dev Ne_sQTL]$ diff -y SRR_freezeSORT_bam.txt convertedbams.txt 
...
SRR5125340.sra.bam					      |	SRR2164650.sra.bam
SRR5125397.sra.bam					      |	SRR2167857.sra.bam
SRR5125400.sra.bam					      |	SRR2170217.sra.bam
SRR5125580.sra.bam					      |	SRR2170696.sra.bam
SRR5125595.sra.bam					      |	SRR2170864.sra.bam
...
```
The files on the left are those that I have not already converted to `bam` but are required by the analysis freeze.

```
SRR5125340.sra
SRR5125397.sra
SRR5125400.sra
SRR5125580.sra
SRR5125595.sra
```


```
#!/bin/bash
#SBATCH --job-name=sra2bam
#SBATCH --time=16:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=24
#SBATCH --array=1-5

ml samtools
ml sra-tools

line=`sed "${SLURM_ARRAY_TASK_ID}q;d" to_convert_toBam.txt`

time {
	sam-dump $line | samtools view --threads 23 -bS - > ${line}.bam
}
```

After they're done converting to bam, I have to filter and convert them to junc: `Reaminingto_filter.txt`

### 02/04/2019
#### Mon 04 Feb 2019 11:09:54 AM EST 
I'm going to just delete all lines in `tissue_table.txt` that contain non-whole blood information and see if that does anything.

Okay I just found out `tissue_table.txt` was from the last run-through I did with the other guys. I'm going to download just the whole blood metadata, create `tissue_table.txt` from that and use that for `sraNameChangeSort.R`. Okay it works. Rajiv sent me the analysis freeze of all the samples that were used for the eQTL by GTEx. I just want to note that I have 454 samples but according to dbGaP, there are actually 465 samples. So I don't know what to do about that. Hopefully they're not all necessary.

Research is so messy. Using `comm -3 smplessorted.txt havesorted.txt` to compare the SRR files I **have** to the SRR metadata, I apparently do not have these 11 files:
```
	SRR5125340
	SRR5125341
	SRR5125397
	SRR5125398
	SRR5125400
	SRR5125401
	SRR5125580
	SRR5125581
	SRR5125595
	SRR5125596
	SRR607214
```

So now I have to check to see if they were in the freeze.

I don't know how to see if any of these SRRs that I didn't download are in the analysis freeze or not. I have that annotations file that I would have to use.

```
IFS=$'\n'       # make newlines the only separator
for sample in $(cat remainingDL.txt)
do
    ./prefetch -X 50G --ascp-path '/software/apps/aspera/3.7.2.354/bin/ascp|/software/apps/aspera/3.7.2.354/etc/asperaweb_id_dsa.openssh' $sample
done
```

I'm just going to ask this question online re: how to use the annotation file.

I also have to play with this command a bit to get the information I want:

`cat GTEx_v7_Annotations_SampleAttributesDS.txt | sed -e1,1d | awk -F '\t' '{if ($17 == "RNASEQ") print $1}' > gtex_v7_analysis_freeze.txt`

#### Mon 04 Feb 2019 02:11:31 PM EST 
Just met with Rajiv. I should go ahead and download those above files, and basically I'll just figure out for myself how to use the analysis freeze samples.

#### 6:31 PM 2/4/2019
There is no column in the analysis freeze for SRR number, which means that I have to get the GTEx IDs for each of the SRA's, compare them to GTEx IDs as part of the freeze. I don't know.

For the dbGaP metadata, Column 25 is the GTEx sample ID. Column 21 is SRR. 

`while read line; do cut -f1 -d '.'; done < sraDL.txt > SRRs.txt`

Above is how I got all of the SRR ID's of the files I have now. Now I need to get the full GTEx ID's for each of these guys using the whole blood metadata, and lastly see which of those sample IDs were **not** excluded from the eQTL analysis, then convert the GTEx samples that were not excluded back into SRR IDs, then move the corresponding SRRs to a seperate folder and do LC again. Fortunately, except for those 11 SRAs I forgot to download, I won't need to do much conversions.

#### 8:22 PM 2/4/2019

`awk '{ print $21, "\t",  $26 }' WholeBlood.txt > SRRGTEX.txt`

I got the whole blood and SRR metadata.

#### 8:47 PM 2/4/2019
~~I'm stuck. I have `awk 'FNR==NR {a[$1]=$0; next}; $1 in a {print a[$1]}' SRRGTEX.txt SRRs.sorted`. `SRRGTEX` is the metadata, SRRs.sorted is what I have. I'm trying to get the GTEX ID's of the samples I have, match them to the GTEx~~

I'm lost. I have the GTEx freeze samples, and I used that line that Rajiv gave me to get all of the samples used in eQTL analysis, but `wc -l` reveals that that's 11688 entries, which is all of the samples sequenced. All I want is the SRR's of the samples used in eQTL analysis. I can do this by finding the GTEx ID for all samples that are not excluded according to the annotation file, then matching those IDs to the corresponding SRRs, and then using those SRRs (which are like 367 or something) to determine which of the SRR files need to be converted/LC'd. Gotta meet with Rajiv tomorrow.

### 02/03/2019
#### Sun 03 Feb 2019 08:30:24 AM EST 
~~`filter_bam` is done running and quickcheck returned no errors.~~ `bam2junc` is done running.

`interact -p shared -m 12G -t 180:0` for intron clustering.

Done with intron clustering. Maybe do it on the dev node next time. Nothing bad happened but it could have.

Submitted `QTLtools_filter.sh`.

#### Sun 03 Feb 2019 10:36:52 PM EST 
Just ran `sraNameChangeSort.R`. Keep getting outputs of tissues other than just whole blood: amygdala, Lung, skin (sun-exposed), caudate, spinal cord. Not sure why this is; should only be whole blood. I will investiage this tomorrow.

### 02/02/2019
#### Sat 02 Feb 2019 08:26:52 AM EST 
the `*.filt` files are done processing. The slurm files indicate no error but running `samtools quickcheck -v *filt` says all the files have no header like that's a bad thing. I'm going to keep going with the process. 

~~I have to call `bam2junccall.sh` for the next step and I'm going to keep it as it is in terms of loading the modules as opposed to running the executables in the working directory. The reason for this is that I don't want to alter the LeafCutter script too much and make bugs, but also I don't know how LC handles these programs downstream. I don't want to slow down MARCC though so in `bam2junccall.sh` I made it so that only up to 25 jobs can run concurrently.~~ Nah you know what I'm just going to change LC to call the executable in its full path. `/scratch/groups/rmccoy22/Ne_sQTL/sra/samtools-1.9/samtools`. I'm still going to have to load python 2.7-anaconda though. Let's see if it works.

The script ran fine, but now, all of the jobs gave me this warning: 
```
Converting SRR1345412.sra.bam.filt to junc
[E::sam_parse1] missing SAM header
[W::sam_read1] Parse error at line 1
[main_samview] truncated file.
```
I guess something in the filtering step broke the files? because the `*.bam` were fine. 

So I think it has to do something with the `*filt` files no longer being in bam format anymore? I can just view them with `head`, whereas the bam files I have to do `samtools view <bam> | head`. Don't know what to do.

I figured out what the problem is. `samtools view -b` flag outputs as bed. I needed that when filtering those bad boys.

#### Sat 02 Feb 2019 06:27:33 PM EST 
The filter is done. Converting to junc now.

### 02/01/2019
#### Fri 01 Feb 2019 08:24:27 AM EST 
The conversion is complete. Need to validate bam files now.

`samtools quickcheck -v *bam`

```
SRR1333462.sra.bam was missing EOF block when one should be present.
SRR1333462.sra.bam
SRR607018.sra.bam had no targets in header.
SRR607018.sra.bam
```

Only 2? That's impressive.

`./sam-dump SRR1333462.sra | ./samtools view --threads 23 -bS - > "SRR1333462.sra.bam"`
`./sam-dump SRR607018.sra | ./samtools view --threads 23 -bS - > "SRR607018.sra.bam"`

I just ran these two commands on the dev node. I'll check back in a few hours to make sure they went through.

#### Fri 01 Feb 2019 02:59:06 PM EST 
`SRR1333462.sra.bam` is not yet done, but all other bams are complete and valid according to `samtools quickcheck`.

#### Fri 01 Feb 2019 06:55:55 PM EST 
`SRR1333462.sra.bam` is still not done, but I think it's almost done because the original is like 6.4Gb and the bam is like a little more the 6.7 right now.

#### Fri 01 Feb 2019 09:09:55 PM EST 
Done. Going to filter the bams now. Changing `filter_bam.sh` on MARCC to not all load `samtools` as a module because MARCC can't handle it. Gonna call the executable instead, like I did for the conversion. Also setting `--threads 24`.

### 01/31/2019
#### Thu 31 Jan 2019 11:56:16 AM EST 
~~Very weird turn of events: I saw that by the time everything will be done, there will be *454* completed bams, but I forgot that the original number, *447* was the one without the sra's above 20Gb, not 454, so I thought I was missing some sra's so I extracted all of the SRR runs from~~ 

I just ran around in a dumbass circle. **454** is the complete number of SRAs that we are doing. If you have any doubt, call `sacct --starttime 2019-01-28`. The jobs are finishing up. There are few of them left. Going to run `picard-tools` to validate them.

### 01/30/2019
#### 9:18 AM 1/30/2019
One of our jobs, 32327199_80, is taking 20+ hours to process, which is unusual because line 80 correpsonds to SRR1333462.sra in sralist.txt, so this means something strange is happening here. I only allotted 24 hours to each job so if it hits the wall, I might have to do it again manually.

I reduced the time to 16 hours and job 80 hit the wall. I'm going to individually convert it on the dev node.

```
time {
        ./sam-dump SRR1333462.sra | samtools view --threads 23 -bS - > SRR1333462.sra.bam
}
```

#### Wed 30 Jan 2019 03:47:30 PM EST 
We got a ton of failed jobs that I have to start from the beginning. We got two different types of error messages and I ascertained which jobs they were using `grep "error" slurm-32327199_* | awk -F'[_.]' '{print $2}' > errorfailure.txt` and `grep "Home" slurm-32327199_* | awk -F'[_.]' '{print $2}' > homedirfailure.txt`. I determined they have no lines in common and concatenated them. I'm going to figure out how to make each line's number correspond to a line in another file and then replace them with the SRA file name. I'm also going to change `sra2bam.sh` not to use the `sra-tools` module and instead utilize a bin. 

I replaced all of the numbers I concatenated into `failedjobs.txt` with the corresponding line in `sralist.txt` into a file `failedsras.txt`.
`awk 'NR==FNR {a[FNR]=$0;next} {printf "%s\t%s\n", a[$1], $2}' sralist.txt failedjobs.txt > failedsras.txt`

I've saved a new file called `failed_sra2bam.sh`.

#### Wed 30 Jan 2019 10:02:26 PM EST 
I just tried to run this job because MARCC finally started working but the script would create files with filenames like `SRR660378.sra?.bam`. Using `cat -A` revealed that the names were actually like `SRR660378.sra^I$`, so I did this: `cat -A failedsras.txt | sed 's/...$//' > failedsrascorrected.txt` to fix it. Now I have to update the script `failed_sra2bam.sh`. Let's see if it works. It works.

### 01/29/2019
#### Tue 29 Jan 2019 10:27:19 AM EST 
I've changed `sra2bam.sh` last night to both use 24 cores (`samtools` has multi-threaded capabilities) and 24 hours. This might be overkill, so I've set up a `time` wrapper to give us an idea of how long these jobs take. I'm guessing that both using 24 cores and 24 hours is too much, I could probably get away with just 16 hours, since the 50-or-so jobs that have finished so far have taken max 5 hours. That would help with the punishingly long wait in the priority queue. But I have time and I don't know how long these files will actually take to convert, just like I didn't know 20Gb was not enough to download all of the SRAs, so I'll just keep it this way and remember to tweak the script once I find out how long it really takes.

Guy named Adam emailed me about something, needs *MY* help.

Have to 1) optimize LeafCutter for use on a bajillion files and 2) figure out the deal with "Picard-tools"

Found out how to use [picard-tools](https://broadinstitute.github.io/picard/command-line-overview.html#ValidateSamFile):

```
java -jar picard.jar ValidateSamFile \
      I=input.bam \
      MODE=SUMMARY
```

re: LeafCutter steps, I already addressed this earlier in the log. I should just call an `interact` with a good amount of memory, maybe like 12 Gb or something.

### 01/28/2019
#### 9:08 AM 1/28/2019
We have the bam files to work with. Going down the list. 

I have to rmember to change some of the array ranges on these batch scripts.

A lot of the files were giving me the following error:
```
Filtering SRR1376294.sra.bam
[W::bam_hdr_read] EOF marker is absent. The input is probably truncated
[E::bgzf_read] Read block operation failed with error 4 after 0 of 4 bytes
[main_samview] truncated file.
```
Not sure what that is so I'm going to keep moving on.

Turns out that meant that the files were incomplete, so I'm re-converting sra -> bam. In doing so, three files encountered errors:

25th sra in sralist.txt...
```
slurmstepd: error: Munge decode failed: Expired credential
slurmstepd: error: Verifying authentication credential: Expired credential
slurmstepd: error: *** JOB 32321214 ON compute0221 CANCELLED AT 2019-01-28T10:43:09 DUE TO JOB REQUEUE ***
```

292nd...
```
/software/lmod/lmod/init/bash: line 70: 26150 Killed                  singularity exec /software/apps/sra-tools/2.9.2/sra-tools-2.9.2.simg sam-dump "$@"
[W::sam_read1] Parse error at line 44952410
[main_samview] truncated file.
slurmstepd: error: Detected 1 oom-kill event(s) in step 32321481.batch cgroup. Some of your processes may have been killed by the cgroup out-of-memory handler.
```

48th...
```
/software/lmod/lmod/init/bash: line 70:  8241 Killed                  singularity exec /software/apps/sra-tools/2.9.2/sra-tools-2.9.2.simg sam-dump "$@"
[W::sam_read1] Parse error at line 13112870
[main_samview] truncated file.
slurmstepd: error: Detected 1 oom-kill event(s) in step 32321237.batch cgroup. Some of your processes may have been killed by the cgroup out-of-memory handler.
```

I tried firing them off individually on the command line. Will report back. **Definitely need to develop a contingency plan if any of the files ever run into errors like this.**

```
[aseyedi2@jhu.edu@compute0230 sra]$ sam-dump SRR603421.sra | samtools view -bS - > SRR603421.sra.bam
[E::sam_parse1] SEQ and QUAL are of different length
[W::sam_read1] Parse error at line 5402763
[main_samview] truncated file.
```

```
[aseyedi2@jhu.edu@compute0230 sra]$ sam-dump SRR1315517.sra | samtools view -bS - > SRR1315517.sra.bam
[W::sam_read1] Parse error at line 13112870
[main_samview] truncated file.
```

I'm still having issues. MARCC crashed so I couldn't finish the last one's conversion (25th).

I'm going to `sha1sum` and re-DL the files to see if the files were corrupted during download

##### 2:47 PM 1/28/2019
I'm still wating for all of the SRAs to finish `sha1sum`, but here are the three of interest:

```
42253bad8a7e7e38b3a3eb5c795e6c8c5bf4688c  SRR1090907.sra
3eba112c500ee99b2d42f600c0292d6da5ab58bf  SRR603421.sra
70461faa76e66bd6d1cb8097abdf36d5d6b6146b  SRR1315517.sra
```
Above are from when I downloaded them the other week. I'm downloading them again to make sure they're of the same quality

```
42253bad8a7e7e38b3a3eb5c795e6c8c5bf4688c  SRR1090907.sra
3eba112c500ee99b2d42f600c0292d6da5ab58bf  SRR603421.sra
70461faa76e66bd6d1cb8097abdf36d5d6b6146b  SRR1315517.sra
```
They're the same, so either MARCC has a problem or dbGaP/GTEx has a problem. I'm going to try converting `SRR1090907.sra` again since I think that one could have worked.


Added several more cores to `sra2bam.sh` because of the multi-threaded capabilities. Not that it matters now.

I got this error on jobs 164 and 165:
```
slurmstepd: error: Munge decode failed: Expired credential
slurmstepd: error: Verifying authentication credential: Expired credential
```
SRR1398728.sra and SRR1399433.sra. 

Big problem. Some of the whole bloods were not downloaded because I set the size too small (20Gb is too small??) gonna try again.

#### 11:10 PM 1/28/2019

Okay, so really there's 454 if you count the big boys. I updated the appropriate scripts to reflect that, and redid the `ls *sra >> sralist.txt` to contain all SRAs.

Submitted the massive job. I also am making a new `sha1sum` of all of the SRAs for whole blood so I'll have to remember that.


### 01/25/2019
#### Fri 25 Jan 2019 10:27:47 AM EST 
I'm going to run the code again to demonstrate to Rajiv exactly what the problem is. 

ok nvm i guess there is supposed to be 2 of every tissue type.

Alright I just talked to Rajiv and we're just going to make this as easy as possible and not expend so much energy to make it harder for ourselves. 

I'm not even going to include brain samples in our analysis right now because dbGaP isn't working.

There are **447** whole blood samples we are working with.

I'm doing the whole thing again but for the whole blood samples.

```
IFS=$'\n'       # make newlines the only separator
for sample in $(cat WholeBloodSRR.txt)
do
    ./prefetch -X 50G --ascp-path '/software/apps/aspera/3.7.2.354/bin/ascp|/software/apps/aspera/3.7.2.354/etc/asperaweb_id_dsa.openssh' $sample
done
```

#### 10:27 PM 1/28/2019
Tons of errors, and incomplete files. I feel like I have to scrap all the bam files and try again. I cancelled the remaining jobs just now.

### 01/23/2019
#### Wed 23 Jan 2019 03:24:59 PM EST 
Now I have to concatenate stuff.

There's a problem. `tissue_table.txt` that is generated by `sraTissueExtract.R` looks wrong. `tissue_table.txt` has duplicates in it and `sraNameChangeSort.R` spits out every type of tissue instead of just one kind. This sucks, especially since it worked before I switched to using the codes.

### 01/22/2019
#### Tue 22 Jan 2019 10:11:07 AM EST 
Just finished indexing the big boy. 

So much tedium. I have to fix the concatenated phenotype files because the header occurs multiple times throughout the file. I also have to figure out how to call QTLtools on every tissue-phenotype file while including also the covariate. I also just found out that the GTEx consortium is generating sQTLs using LeafCutter and it's undoubtably going to be better than everything I've done here for the past like... 4 months lol. It's fine I'm not upset I just think it's funny. Life. So it makes me wonder how important it is to actually pour so much energy into having a functioning pipeline? like it really doesn't matter, right? since I'm going to only test this out on just one tissue type, which we've decided will be whole blood.

Whatever I'm almost done.

```
ml sra-tools
IFS=$'\n'       # make newlines the only separator
for sample in $(cat WholeBloodSRR.txt)
do
    ./prefetch -X 20G --ascp-path '/software/apps/aspera/3.7.2.354/bin/ascp|/software/apps/aspera/3.7.2.354/etc/asperaweb_id_dsa.openssh' $sample
done
```

That's how I'm downloading all of the whole blood samples

```
for tissue in $(cat ./phenpaths.txt)    
do
    awk '{ if($0 != header) { print; } if(header == "") { header=$0; } }' $tissue | sed 's/^[0-9][0-9]*[ \t]*//' > $tissue
done
```

Now I have to rename everything to some sort of tissue code.


Here they are
```
Tissue	Key
Adipose - Subcutaneous	ADPSBQ
Adipose - Visceral (Omentum)	ADPVSC
Adrenal Gland	ADRNLG
Artery - Aorta	ARTAORT
Artery - Coronary	ARTACRN
Artery - Tibial	ARTTBL
Bladder	BLDDER
Brain - Amygdala	BRNAMY
Brain - Anterior cingulate cortex (BA24)	BRNACC
Brain - Caudate (basal ganglia)	BRNCDT
Brain - Cerebellar Hemisphere	BRNCHB
Brain - Cerebellum	BRNCHA
Brain - Cortex	BRNCTXA
Brain - Frontal Cortex (BA9)	BRNCTXB
Brain - Hippocampus	BRNHPP
Brain - Hypothalamus	BRNHPT
Brain - Nucleus accumbens (basal ganglia)	BRNNCC
Brain - Putamen (basal ganglia)	BRNPTM
Brain - Spinal cord (cervical c-1)	BRNSPC
Brain - Substantia nigra	BRNSNG
Breast - Mammary Tissue	BREAST
Cells - EBV-transformed lymphocytes	LCL
Cells - Transformed fibroblasts	FIBRLBLS
Cervix - Ectocervix	CVXECT
Cervix - Endocervix	CVSEND
Colon - Sigmoid	CLNSGM
Colon - Transverse	CLNTRN
Esophagus - Gastroesophageal Junction	ESPGEJ
Esophagus - Mucosa	ESPMCS
Esophagus - Muscularis	ESPMSL
Fallopian Tube	FLLPNT
Heart - Atrial Appendage	HRTAA
Heart - Left Ventricle	HRTLV
Kidney - Cortex	KDNCTX
Liver	LIVER
Lung	LUNG
Minor Salivary Gland	SLVRYG
Muscle - Skeletal	MSCLSK
Nerve - Tibial	NERVET
Ovary	OVARY
Pancreas	PNCREAS
Pituitary	PTTARY
Prostate	PRSTTE
Skin - Not Sun Exposed (Suprapubic)	SKINNS
Skin - Sun Exposed (Lower leg)	SKINS
Small Intestine - Terminal Ileum	SNTTRM
Spleen	SPLEEN
Stomach	STMACH
Testis	TESTIS
Thyroid	THYROID
Uterus	UTERUS
Vagina	VAGINA
Whole Blood	WHLBLD

```

They're in `scratch` for now. I screwed something up but this presents with a good opportunity to start over. I'm going to start from the post-LC step.


What happened today? I started over from an earlier step. I'm now waiting on an answer on how to write an R script that replaces all full tissue names with the code. That would be great and fix all of my problems. Then I would just go about the rest of the thing and I think I would be done.

I did it. I changed the phenotype files to the abbreviation code. Now I have to do the same with covariates. I feel like Unix is like Latin.


### 01/21/2019
#### Mon 21 Jan 2019 10:26:30 AM EST 
It's the King's day, as my high school wrestling coach once said. Here's a quote in his memory:

```
“When machines and computers, profit motives and property rights are considered more important than people, the giant triplets of racism, extreme materialism and militarism are incapable of being conquered.”

—  “Revolution of Values,” 1967
```

I'm finally back at using QTLtools.

Okay, MARCC crashed and isn't working. Don't know how to fix it. 

Works again. Okay, I'm about to run QTLtools **EXCEPT** I have no idea where I left the GTEx VCF file so I have to copy & filter & index it again like I did in an earlier step, which will take a while. But I'm almost done.

### Mon 21 Jan 2019 04:43:43 PM EST 
I'm running `bcftools` and just waiting on it to finish.


### 01/18/2019
#### Fri 18 Jan 2019 11:06:53 AM EST 
`for i in {1..48}; do line=``sed "${i}q;d" tissuenames.txt``; echo "Concatenating $line..."; for q in {1..22}; do echo "Chr $q..."; cat "${q}_${line}.txt" >> [1-22]${line}.txt; done; done`

`cat tissue_table.txt | cut -f3 | awk '{if(NR>1)print}'`

Turns out that the file names for the covariates are formatted differently than the tissue names that come from the SRA metadata. gonna fix now. 

Works.

### 01/17/2019
#### Thu 17 Jan 2019 03:13:57 PM EST 
Turns out the fine folk over at GTEx consortium are gonna beat us to the punch when it comes to generated sQTL calls in v8. It's okay. Rajiv helped me out with `sraNameChangeSort.R`. Turns out I don't have to do anything with the covariates. I guess I forgot. 

#### Thu 17 Jan 2019 05:07:30 PM EST 
It's getting to be that time of the day where MARCC starts crashing. Here is the line I was trying to hit: `for i in {1..48}; do line=``sed "${i}q;d" tissuenames.txt``; echo "Concatenating $line..."; for q in {1..22}; do echo "Chr $q..."; cat "${q}_${line}.txt" >> [1-22]${line}.txt; done; done`

### 01/16/2019
#### Wed 16 Jan 2019 12:54:55 PM EST 
Our massive VCF was indexed overnight and just now I deleted it thinking that it wasn't fully indexed. Turns out the file doesn't even appear until it's done indexing. I'm going to do a `sha1sum` on it so I have a record of exactly what it's supposed to look like. Done.

#### Wed 16 Jan 2019 02:50:34 PM EST 

```
Warning message:
In fread(args[1]) :
  Stopped early on line 2. Expected 1 fields but found 1. Consider fill=TRUE and comment.char=. First discarded non-empty line: <<>>
Error in `[.data.table`(NE, , .SD, .SDcols = c(keep, intersect(names(NE),  : 
  Some items of .SDcols are not column names (or are NA)
Calls: lapply -> FUN -> [ -> [.data.table
```
I'm getting this error when running `sraNameChangeSort.R` and I don't know why. I should talk to Rajiv about it.

Update: I'm just going to talk to Rajiv about it and take it easy for today. I don't want to screw anything up. 


### 01/15/2019
#### Tue 15 Jan 2019 02:08:25 PM EST 
I finished making `sraNameChange&Sort.R`. Hopefully it works. Rajiv is now working on getting access to the CRAM files hosted on Google/Amazon but I think we're going to test the pipeline again after this with SRAs of just one tissue. Now I have to, very painfully, do the same thing again, but concoct a script that replaces the *subject IDs* that are listed as column names in the GTEX covariate files with the *sample IDs*, which are longer and denote the tissue and subject, instead of just the subject. I was thinking about extracting the GTEX subject IDs and replacing it with the corresponding sample ID **IF** both the tissue type and subject ID match. This is going to be yet another struggle.

To do this, I'm going to have to reconfigure `sraTissueExtract.R` to also print out `submitted_subject_id`. Okay cool. 

#### Tue 15 Jan 2019 03:49:57 PM EST 
I'm having trouble figuring out if, when calling a script outside of the working directory, one has to provide the path for the thing as;kdjf;lksajdf;lkkjs nevermind

#### Tue 15 Jan 2019 04:16:50 PM EST 
I made an oopsie daisy and have to start again.

Okay, I give up. I think MARCC is broken. I tried to retrace my steps to figure out exactly where I was and what I had to do next and I just ended up in prison. I kind of want to start over again from the beginning whenever I could get this to work.

### 01/14/2019
#### Mon 14 Jan 2019 10:54:43 AM EST 
~~After a couple of days of procrastinating, I am ready to face this monster. This demon. First, last night or maybe two nights ago I realized that I should probably rename these guys **before** running `prepare_phenotype.py` or whatever it's called. This is because, when it generates the covariates, it makes the Run ID the colname. I could also change the colnames of the covariate file in place, which is what I think I will do.~~

~~I'm going to try to change everything after running `prepare_phenotype.py`. This is going to be hard. I'm going to download these files so I can develop the R script locally.~~

~~Okay. I'm pretty confused and am not sure what to do. The GTEx supplied covariates are all in subject id format i.e. GTEX-.... (with an optional 5th char). I have a tissue table which has sra run IDs, corresponding tissue type, and *full* GTEx sample ID i.e. GTEX-....-....-..-.... (or something like that). Not sure how to proceed here. Do I want to change the LeafCutter PCs file to match the GTEx supplied one? My puny brain which was intentionally designed for being dumb is having a meltdown.~~

I'm confusing the hell out of myself. I think I should just roll back a step, to before `prepare_phenotype.py`.

Starting over with `python ../../../leafcutter/clustering/leafcutter_cluster.py -j juncfiles.txt -r intronclustering/ -m 50 -o NE_sQTL -l 500000`. 

Okay, I reread the notes I took. I guess I forgot that I was doing this (this being... renaming column headers and segregating based on tissue) for only the phenotype files. But I also have to do this sort of thing for the covariate files. I guess I could acheive this by simply renaming the colnames for the LC-generated PCs and then *renaming the GTEx-supplied covariate colnames* to the sample id of the tissue/subject complex.

This is so complicated. I have no idea how to make this more simple. I do not like how complicated it is. 

I'm just going to go ahead with doing this rename and rewrite thing but only for the phenotype files right now. I will later focus on the PCs.

#### Mon 14 Jan 2019 04:51:25 PM EST 
[Some guy on stack overflow helped me figure it out.](https://stackoverflow.com/questions/54189236/how-do-i-split-a-table-into-several-new-tables-based-on-whether-the-column-heade/54189469#54189469)

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

