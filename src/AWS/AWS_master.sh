##########################################################################################################################################
#  ███╗   ██╗███████╗ █████╗ ███╗   ██╗██████╗ ███████╗██████╗ ████████╗██╗  ██╗ █████╗ ██╗         ███████╗ ██████╗ ████████╗██╗     
#  ████╗  ██║██╔════╝██╔══██╗████╗  ██║██╔══██╗██╔════╝██╔══██╗╚══██╔══╝██║  ██║██╔══██╗██║         ██╔════╝██╔═══██╗╚══██╔══╝██║     
#  ██╔██╗ ██║█████╗  ███████║██╔██╗ ██║██║  ██║█████╗  ██████╔╝   ██║   ███████║███████║██║         ███████╗██║   ██║   ██║   ██║     
#  ██║╚██╗██║██╔══╝  ██╔══██║██║╚██╗██║██║  ██║██╔══╝  ██╔══██╗   ██║   ██╔══██║██╔══██║██║         ╚════██║██║▄▄ ██║   ██║   ██║     
#  ██║ ╚████║███████╗██║  ██║██║ ╚████║██████╔╝███████╗██║  ██║   ██║   ██║  ██║██║  ██║███████╗    ███████║╚██████╔╝   ██║   ███████╗
#  ╚═╝  ╚═══╝╚══════╝╚═╝  ╚═╝╚═╝  ╚═══╝╚═════╝ ╚══════╝╚═╝  ╚═╝   ╚═╝   ╚═╝  ╚═╝╚═╝  ╚═╝╚══════╝    ╚══════╝ ╚══▀▀═╝    ╚═╝   ╚══════╝
##########################################################################################################################################

# adjust for the interactive session time length you might need here

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

# this project's scripts dir
scripts=$(echo /home-1/aseyedi2@jhu.edu/work/aseyedi2/neand_sQTL/src/primary/)
# data dir
data=$(echo /home-1/aseyedi2@jhu.edu/work/aseyedi2/neand_sQTL/data/)
# ncbi/files/
ncbiFiles=$(echo /scratch/groups/rmccoy22/Ne_sQTL/files/)
# IF YOU ALREADY HAVE NON-BIALLELIC INDEXED VCF
VCF=$(echo /scratch/groups/rmccoy22/Ne_sQTL/files/GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz)
# input directory with junc files here
junc=$(echo '/scratch/groups/rmccoy22/Ne_sQTL/sra/sqtl_junc')
# leafcutter directory here
leafCutter=$(echo /scratch/groups/rmccoy22/aseyedi2/leafcutter)

mkdir intronclustering/

sbatch --wait ${scripts}/../AWS/junc_cluster.sh

cd intronclustering/

sbatch --wait ${scripts}/../AWS/prepare_phen_table.sh

## Step 4 - VCF Preparation (optional, see doc for details)
################################################

## Step 5 - QTLtools Preparation
################################################
# prepare files for QTLtools
ls *qqnorm*.gz | sort -V > leafcutterphenotypes.txt 
# important: render these files compatible with QTLtools
echo "Making phenotype files QTLtools compatible..."
sbatch --wait ${scripts}/sh/QTLtools-Filter.sh
ls *.qtltools | sort -V > qtltools-input.txt
# generate the corresponding tbi files
rm Ne*tbi
for i in {1..22}; do tabix -p bed Ne-sQTL_perind.counts.gz.qqnorm_chr${i}.gz.qtltools; echo "Bedding chromosome $i"; done

cp ${data}/01-22-2019/GTExTissueKey.csv $PWD
# get the tissue sites for each corresonding sra file
Rscript ${scripts}/R/sraTissueExtract.R ${data}/Metadata/SraRunTable.txt GTExTissueKey.csv

# submit each LF phenotype file to sraNameChangeSort as command line variable as well as tissue_table.txt
for phen in *qqnorm*.gz.qtltools; do Rscript ${scripts}/R/sraNameChangeSort.R $phen tissue_table.txt ; done
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

cp $data/Metadata/GTExCovKey.csv $PWD

# Moves covariates to corresponding directory
for line in $(cat GTExCovKey.csv)
do
   full=$(echo $line | awk -F',' '{print $1}')
   abb=$(echo $line | awk -F',' '{print $2}')
   if grep "$abb" tissuesused.txt; then
      cp GTEx_Analysis_v7_eQTL_covariates/$full.v7.covariates.txt $abb
      Rscript ${scripts}/R/mergePCs.R Ne-sQTL_perind.counts.gz.PCs $abb/$full.v7.covariates.txt tissuetable/tissue_table.txt
      mv $full.v7.covariates_output.txt $abb
   fi
done

## Step 4 - Mapping sQTLs using QTLtools
################################################
for line in $(cat GTExCovKey.csv); do
   full=$(echo $line | awk -F',' '{print $1}')
   abb=$(echo $line | awk -F',' '{print $2}')
   if grep "$abb" tissuesused.txt; then

      sbatch --export=abb=$abb,full=$full ${scripts}/sh/QTLTools-Loop.sh
      
   fi
done

sbatch ${scripts}/sh/QQViz.sh