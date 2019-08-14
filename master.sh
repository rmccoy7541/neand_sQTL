#!/bin/bash
##########################################################################################################################################
#  ███╗   ██╗███████╗ █████╗ ███╗   ██╗██████╗ ███████╗██████╗ ████████╗██╗  ██╗ █████╗ ██╗         ███████╗ ██████╗ ████████╗██╗        #
#  ████╗  ██║██╔════╝██╔══██╗████╗  ██║██╔══██╗██╔════╝██╔══██╗╚══██╔══╝██║  ██║██╔══██╗██║         ██╔════╝██╔═══██╗╚══██╔══╝██║        #
#  ██╔██╗ ██║█████╗  ███████║██╔██╗ ██║██║  ██║█████╗  ██████╔╝   ██║   ███████║███████║██║         ███████╗██║   ██║   ██║   ██║        #
#  ██║╚██╗██║██╔══╝  ██╔══██║██║╚██╗██║██║  ██║██╔══╝  ██╔══██╗   ██║   ██╔══██║██╔══██║██║         ╚════██║██║▄▄ ██║   ██║   ██║        #
#  ██║ ╚████║███████╗██║  ██║██║ ╚████║██████╔╝███████╗██║  ██║   ██║   ██║  ██║██║  ██║███████╗    ███████║╚██████╔╝   ██║   ███████╗   #
#  ╚═╝  ╚═══╝╚══════╝╚═╝  ╚═╝╚═╝  ╚═══╝╚═════╝ ╚══════╝╚═╝  ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝    ╚══════╝ ╚══▀▀═╝    ╚═╝   ╚══════╝   #
##########################################################################################################################################
# Fill out the variables below with the appropriate full paths for each of the corresponding directories.                                #
# It is also assumed that the GTEx VCF file for the whole genome has already been downloaded (and processed)                        	 #
# 	(see Documentation for details)																										 #
##########################################################################################################################################

# load modules
ml bedtools
ml samtools
ml sra-tools
ml python/2.7-anaconda
ml bcftools
ml htslib
ml R
ml gcc
ml

#### Fill the Directories appropriately here
# this project's scripts dir
scripts=$(echo {path/to}/neanderthal-sqtl/src/primary/)
# data dir
data=$(echo {path/to}/neanderthal-sqtl/data/)
# IF YOU ALREADY HAVE NON-BIALLELIC INDEXED VCF
VCF=$(echo {path/to}/GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz)
# leafcutter directory here
leafCutter=$(echo {path/to}/leafCutter/)
# sprime
sprime=$(echo ${data}/SPRIME/sprime_calls.txt)
# base working dir
basewd=$(echo )
# back up dir - there is where the results will go
backupdir=$(echo )

cd $basewd

echo "Clustering introns..."
sbatch --wait --export=LC=$(echo $leafCutter) ${scripts}/sh/01_junc_cluster.sh

cd intronclustering/

echo "Preparing phenotype table..."
sbatch --wait ${scripts}/sh/03_prepare_phen_table.sh $leafCutter

## Step 5 - QTLtools Preparation
################################################
# prepare files for QTLtools
ls *qqnorm* > leafcutterphenotypes.txt 
# important: render these files compatible with QTLtools
echo "Making phenotype files QTLtools compatible..."
sbatch --wait ${scripts}/sh/04_QTLtools-Filter.sh
ls *.qtltools >> qtltools-input.txt

# generate the corresponding tbi files
for i in {1..22}; do tabix -p bed Ne-sQTL_perind.counts.gz.qqnorm_chr${i}.gz.qtltools; echo "Bedding chromosome $i"; done

cp ${data}/GTExTissueKey.csv $PWD
# get the tissue sites for each corresonding sra file
Rscript ${scripts}/R/05_sraTissueExtract.R ${data}/SraRunTable.txt GTExTissueKey.csv

# submit each LF phenotype file to sraNameChangeSort as command line variable as well as tissue_table.txt
for phen in *qqnorm*.gz.qtltools; do Rscript ${scripts}/R/06_sraNameChangeSort.R $phen tissue_table.txt ; done
cat tissue_table.txt | cut -f3 | awk '{if(NR>1)print}' |  awk '!seen[$0]++' > tissuenames.txt

mkdir tissuetable/
mv tissue_table.txt tissuetable/
# make directories for each type of tissue
for i in 1_*.txt; do echo $i | cut -d'_' -f 2| cut -d'.' -f 1 | xargs mkdir; done

# save tissue types in a file
for i in 1_*.txt; do echo $i | cut -d'_' -f 2| cut -d'.' -f 1 >> tissuesused.txt; done

# moves each outputted file into its respective tissue folder
for i in *_*.txt; do echo $i | awk -F'[_.]' '{print $2}' | xargs -I '{}' mv $i '{}' ; done

## Concatting the phenotype files
for line in $(cat tissuesused.txt)
do
   head -1 $line/1_$line.txt > $line/$line.phen_fastqtl.bed
   echo "Concatenating $line phenotypes..."
   for file in $(ls $line/*_*.txt)
   do
      echo "Adding $file..."
      cat $file | sed -e1,1d >> $line/$line.phen_fastqtl.bed
   done
done

for line in $(cat tissuesused.txt)
do
   echo "Sorting $line.phen_fastqtl.bed to $line/$line.pheno.bed..."
   bedtools sort -header -i $line/$line.phen_fastqtl.bed > $line/$line.pheno.bed 
   echo "bgzipping $line/$line.pheno.bed..."
   bgzip -f $line/$line.pheno.bed
   #figure out where tabix outputs
   echo "Indexing $line/$line.pheno.bed.gz..."
   tabix -p bed $line/$line.pheno.bed.gz
done

for line in $(cat tissuesused.txt)
do
   mkdir $line/sepfiles
   mv $line/*_${line}.txt $line/sepfiles/
done

# download genotype covariates
wget https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL_covariates.tar.gz

tar -xvf GTEx_Analysis_v7_eQTL_covariates.tar.gz

cp $data/GTExCovKey.csv $PWD

# Moves covariates to corresponding directory
for line in $(cat GTExCovKey.csv)
do
   full=$(echo $line | awk -F',' '{print $1}')
   abb=$(echo $line | awk -F',' '{print $2}')
   if grep "$abb" tissuesused.txt; then
      cp GTEx_Analysis_v7_eQTL_covariates/$full.v7.covariates.txt $abb
      Rscript ${scripts}/R/07_mergePCs.R Ne-sQTL_perind.counts.gz.PCs $abb/$full.v7.covariates.txt tissuetable/tissue_table.txt
      mv $full.v7.covariates_output.txt $abb
   fi
done

## Step 4 - Mapping sQTLs using QTLtools
################################################
numTissues=$(wc -l GTExCovKey.csv)

# Will take at least 3 weeks lol
sbatch --wait \
  --export=scripts=$scripts,data=$data,vcf=$vcf,sprime=$sprime \
  -a 2-$numTissues \
  ${scripts}/sh/08_QTLTools-Loop.sh 

mkdir -p $backupdir
mkdir $backupdir/all_noms
mkdir $backupdir/sqtl_nrich
mv $PWD/*/*permutation* $backupdir
mv $PWD/*/*pheno* $backupdir
mv tissuenames.txt $backupdir

cd $backupdir
cp $sprime .

# concatenates permutation pass results
for TISSUE in ADPSBQ ADPVSC ADRNLG ARTACRN ARTAORT ARTTBL BREAST BRNACC BRNAMY BRNCDT BRNCHA BRNCHB BRNCTXA BRNCTXB BRNHPP BRNHPT BRNNCC BRNPTM BRNSNG BRNSPC CLNSGM CLNTRN ESPGEJ ESPMCS ESPMSL FIBRLBLS HRTAA HRTLV LCL LIVER LUNG MSCLSK NERVET OVARY PNCREAS PRSTTE PTTARY SKINNS SKINS SLVRYG SNTTRM SPLEEN STMACH TESTIS THYROID UTERUS VAGINA WHLBLD
do
  mkdir ${TISSUE}
  mv ${TISSUE}*permutations* ${TISSUE}
  for i in {1..100}
  do
    echo "Catting ${TISSUE}_permutations_chunk_${i}.txt..."
    cat ${TISSUE}/${TISSUE}_permutations_chunk_${i}.txt >> ${TISSUE}_permutations.txt
  done
done

# getting tissue names
for i in $(ls *_permutations.txt | sort -V); do echo $i | cut -d'_' -f 1; done > tissuenames.txt

wget ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Assembly/GRCh37-HG19_Broad_variant/Homo_sapiens_assembly19.fasta

# get allele frequences from VCF
sbatch --export=VCF=$VCF $scripts/sh/13_VariantToTable.sh

# concat the GTEx AF VCF chunks
cat GTExWGS.AF.chr1.txt > GTExWGS.AF.all.txt
for i in {2..22}; do tail -n +2 GTExWGS.AF.chr${i}.txt >> GTExWGS.AF.all.txt ; done
rm GTExWGS.AF.chr*.txt