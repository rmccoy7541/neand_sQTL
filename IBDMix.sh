#!/bin/bash
ml gcc
ml cmake
ml htslib
ml plink
ml bcftools
mkdir /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/vcfChrom
# minimize VCFs while preserving ref/alt allele encoding
for i in {1..22}
do
plink \
  --vcf /scratch/groups/rmccoy22/aseyedi2/NL_sQTL_iso/IBDmix/archaicGenomes/vcfChrom/AltaiNea.hg38_1000g.${i}.LowQualRemoved.vcf \
  --keep-allele-order \
  --real-ref-alleles \
  --recode vcf \
  --allow-extra-chr \
  --out /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/vcfChrom/altai_plink_recode_chr${i} &
done
wait
plink \
  --vcf /scratch/groups/rmccoy22/aseyedi2/NL_sQTL_iso/IBDmix/archaicGenomes/GTExWGSGenotypeMatrixBiallelicOnly_v8.vcf \
  --keep-allele-order \
  --real-ref-alleles \
  --recode vcf \
  --out /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/vcf/gtex_plink_recode
# bgzip and index the lifted over archaic VCFs
for i in {1..22}
do
  (bgzip /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/vcfChrom/altai_plink_recode_chr${i}.vcf ;
  tabix -p vcf /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/vcfChrom/altai_plink_recode_chr${i}.vcf.gz) &
done
wait
# sort the archaic VCFs
for i in {1..22}
do
  (mkdir /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/vcfChrom/tmp_chr${i} ;
  bcftools sort \
    --output-type z \
    --temp-dir /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/vcfChrom/tmp_chr${i} \
    --max-mem 20G \
    --output-file /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/vcfChrom/altai_plink_recode_sorted_chr${i}.vcf.gz \
  /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/vcfChrom/altai_plink_recode_chr${i}.vcf.gz ;
  tabix -p vcf /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/vcfChrom/altai_plink_recode_sorted_chr${i}.vcf.gz) &
done
wait
# combine the per-chromosome archaic VCFs (as a minorit of records still belong to other chromosomes)
bcftools concat \
  --allow-overlaps \
  --output-type z \
  --threads 48 \
  --output /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/vcf/altai_plink_recode_combined_unsorted.vcf.gz \
  /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/vcfChrom/altai_plink_recode_sorted_chr*.vcf.gz
tabix -p vcf /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/vcf/altai_plink_recode_combined_unsorted.vcf.gz
# then split them again
for i in {1..22}
do
  bcftools view \
    -r ${i} /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/vcf/altai_plink_recode_combined_unsorted.vcf.gz \
    > /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/vcfChrom/altai_generate_gt_input_chr${i}.vcf &
done
wait
bgzip -@ 48 /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/vcf/gtex_plink_recode.vcf
tabix -p vcf /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/vcf/gtex_plink_recode.vcf.gz
for i in {1..22}
do
  bcftools view \
    -r ${i} /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/vcf/gtex_plink_recode.vcf.gz \
    > /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/vcfChrom/gtex_generate_gt_input_chr${i}.vcf &
done
wait
for i in {1..22}
do
  /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/IBDmix/build/src/generate_gt \
  --modern /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/vcfChrom/gtex_generate_gt_input_chr${i}.vcf \
  --archaic /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/vcfChrom/altai_generate_gt_input_chr${i}.vcf \
  --output /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/generate_gt_out/altai_gtex_combined_gt_chr${i}.txt &
done
wait
for i in {1..22}
do
  /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/IBDmix/build/src/ibdmix \
  --genotype /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/generate_gt_out/altai_gtex_combined_gt_chr${i}.txt \
  --write-snps \
  --output /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/ibdmix_out/gtex_ibdmix_chr${i}.txt &
done
wait
for i in {1..22}
do
  /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/IBDmix/src/summary.sh \
  50000 \
  4 \
  gtex \
  /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/ibdmix_out/gtex_ibdmix_chr${i}.txt \
  /work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/ibdmix_summary/gtex_ibdmix_chr${i}_50000_lod4.txt &
done
