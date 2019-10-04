#!/bin/bash
#SBATCH --job-name=sprime_combine
#SBATCH --time=4:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=12

### 01-sprime-combine.sh ###

ml bcftools
ml vcftools
ml htslib
ml java

# gtex-sprime directory
wrkdir=/scratch/groups/rmccoy22/aseyedi2/sQTLv8/sprime_test
# 1kg project direcotry (1kg-hs37d5/)
kg_dir=/work-zfs/rmccoy22/resources/1kg-hs37d5/
# original gtex vcf (prefilter) GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz
gtex_vcf=

cd ${wrkdir}

mkdir -p genetic_map
wget http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip
unzip plink.GRCh37.map.zip
mv plink.chr*.GRCh37.map genetic_map/
rm plink.GRCh37.map.zip

# concatenate genetic maps
for i in {1..22}
do 
  cat genetic_map/plink.chr${i}.GRCh37.map | awk '{print "chr"$0}' >> genetic_map/plink.all_autosomes.GRCh37.map
done

# get a list of merged, filtered VCFs to concatenate
for i in {1..22}
do
  echo "${wrkdir}/filtered_vcf/merged_filtered_chr${i}.vcf.gz" >> processed_vcf_list.txt
done

# concatenate VCFs
bcftools concat \
  --file-list processed_vcf_list.txt \
  --naive \
  --threads 12 \
  -O z \
  -o merged_filtered_all_autosomes.vcf.gz

# index
tabix -p vcf merged_filtered_all_autosomes.vcf.gz

touch sprime_v7_combine_complete