#!/bin/bash
#SBATCH --job-name=anno_splice_sites
#SBATCH --time=6:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4

# original gtex vcf (prefilter) GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz
vcf=$(echo )
# path to vcfanno binary
wget https://github.com/brentp/vcfanno/releases/download/v0.3.2/vcfanno_linux64
vcfanno=$(echo vcfanno_linux64)
# path to ensembl vep bin
vep=$(echo )
# vep cache dir
vep_cache=$(echo )

ml htslib
ml bedtools
ml gatk

mkdir spliceai/
cd spliceai/

# use tabix to rapidly extract variants based on overlapping coordinates

split -l 16842 tabix_query.txt tabix_query_

for i in tabix_query_*
do
  while read LINE
  do 
    tabix $vcf ${LINE} >> ${i}.vcf
  done < ${i} &
done

wait

zcat $vcf | head -10000 | grep '^#' > neand_snps_unfiltered.vcf

cat tabix_query_*.vcf >> neand_snps_unfiltered.vcf

bedtools sort -header -i neand_snps_unfiltered.vcf > neand_snps_sorted.vcf

# further filter by requiring matching alleles

gatk SelectVariants \
  --variant neand_snps_sorted.vcf \
  --keep-ids selectvariants_query.list \
  --output neand_snps.vcf

# annotate with SpliceAI and VEP

$vcfanno \
  template.toml \
  neand_snps.vcf > neand_snps_spliceai.vcf

$vep \
  -i neand_snps.vcf \
  --format vcf \
  -o neand_snps_vep.vcf \
  --canonical \
  --port 3337 \
  --symbol \
  --cache \
  --dir_cache $vep_cache \
  --force_overwrite
