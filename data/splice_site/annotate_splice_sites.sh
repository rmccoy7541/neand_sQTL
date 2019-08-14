#!/bin/bash

ml htslib
ml bedtools
ml gatk

cd /work-zfs/rmccoy22/rmccoy22/sqtl/spliceai


# use tabix to rapidly extract variants based on overlapping coordinates

split -l 16842 tabix_query.txt tabix_query_

for i in tabix_query_*
do
  while read LINE
  do 
    tabix /scratch/groups/rmccoy22/Ne_sQTL/files/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz ${LINE} >> ${i}.vcf
  done < ${i} &
done

wait

zcat /scratch/groups/rmccoy22/Ne_sQTL/files/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz | head -10000 | grep '^#' > neand_snps_unfiltered.vcf

cat tabix_query_*.vcf >> neand_snps_unfiltered.vcf

bedtools sort -header -i neand_snps_unfiltered.vcf > neand_snps_sorted.vcf

# further filter by requiring matching alleles

/work-zfs/rmccoy22/progs/gatk-4.1.2.0/gatk \
  SelectVariants \
  --variant neand_snps_sorted.vcf \
  --keep-ids selectvariants_query.list \
  --output neand_snps.vcf

# annotate with SpliceAI and VEP

~/work/progs/vcfanno_linux64 \
  template.toml \
  neand_snps.vcf > neand_snps_spliceai.vcf

~/data/progs/ensembl-vep/vep \
  -i neand_snps.vcf \
  --format vcf \
  -o neand_snps_vep.vcf \
  --canonical \
  --port 3337 \
  --symbol \
  --cache \
  --dir_cache /work-zfs/rmccoy22/rmccoy22/sqtl/spliceai/vep_cache/ \
  --force_overwrite
