#!/bin/bash

ml gcc
ml cmake
ml htslib
ml plink
ml bcftools
ml R

wrkdir="/work-zfs/rmccoy22/rmccoy22/sqtl/ibdmix/"

mkdir -p ${wrkdir}/vcfChrom

# minimize VCFs while preserving ref/alt allele encoding
for i in {1..22}
do
plink \
  --vcf /scratch/groups/rmccoy22/aseyedi2/NL_sQTL_iso/IBDmix/archaicGenomes/vcfChrom/AltaiNea.hg38_1000g.${i}.LowQualRemoved.vcf \
  --keep-allele-order \
  --real-ref-alleles \
  --recode vcf \
  --allow-extra-chr \
  --out ${wrkdir}/vcfChrom/altai_plink_recode_chr${i} &
done

wait

plink \
  --vcf /scratch/groups/rmccoy22/aseyedi2/NL_sQTL_iso/IBDmix/archaicGenomes/GTExWGSGenotypeMatrixBiallelicOnly_v8.vcf \
  --keep-allele-order \
  --real-ref-alleles \
  --recode vcf \
  --out ${wrkdir}/vcf/gtex_plink_recode

# bgzip and index the lifted over archaic VCFs

for i in {1..22}
do
  (bgzip ${wrkdir}/vcfChrom/altai_plink_recode_chr${i}.vcf ;
  tabix -p vcf ${wrkdir}/vcfChrom/altai_plink_recode_chr${i}.vcf.gz) &
done

wait

# sort the archaic VCFs

for i in {1..22}
do
  (mkdir ${wrkdir}/vcfChrom/tmp_chr${i} ;
  bcftools sort \
    --output-type z \
    --temp-dir ${wrkdir}/vcfChrom/tmp_chr${i} \
    --max-mem 20G \
    --output-file ${wrkdir}/vcfChrom/altai_plink_recode_sorted_chr${i}.vcf.gz \
  ${wrkdir}/vcfChrom/altai_plink_recode_chr${i}.vcf.gz ;
  tabix -p vcf ${wrkdir}/vcfChrom/altai_plink_recode_sorted_chr${i}.vcf.gz) &
done

wait

# combine the per-chromosome archaic VCFs (as a minority of records still belong to other chromosomes)

bcftools concat \
  --allow-overlaps \
  --output-type z \
  --threads 48 \
  --output ${wrkdir}/vcf/altai_plink_recode_combined_unsorted.vcf.gz \
  ${wrkdir}/vcfChrom/altai_plink_recode_sorted_chr*.vcf.gz
  
tabix -p vcf ${wrkdir}/vcf/altai_plink_recode_combined_unsorted.vcf.gz 

# then split them again

for i in {1..22}
do
  bcftools view \
    -r ${i} ${wrkdir}/vcf/altai_plink_recode_combined_unsorted.vcf.gz \
    > ${wrkdir}/vcfChrom/altai_generate_gt_input_chr${i}.vcf &
done

wait

# split the GTEx VCF by chromosome

bgzip -@ 48 ${wrkdir}/vcf/gtex_plink_recode.vcf
tabix -p vcf ${wrkdir}/vcf/gtex_plink_recode.vcf.gz

for i in {1..22}
do
  bcftools view \
    -r ${i} ${wrkdir}/vcf/gtex_plink_recode.vcf.gz \
    > ${wrkdir}/vcfChrom/gtex_generate_gt_input_chr${i}.vcf &
done

wait

for i in {1..22}
do
  ~/code/IBDmix/build/src/generate_gt \
  --modern ${wrkdir}/vcfChrom/gtex_generate_gt_input_chr${i}.vcf \
  --archaic ${wrkdir}/vcfChrom/altai_generate_gt_input_chr${i}.vcf \
  --output ${wrkdir}/generate_gt_out/altai_gtex_combined_gt_chr${i}.txt &
done

wait

# modified IBDmix to output single-SNP LOD
# scores for all SNPs that positively contribute to a haplotype
for i in {1..22}
do
  ~/code/IBDmix/build/src/ibdmix \
  -g ${wrkdir}/generate_gt_out/altai_gtex_combined_gt_chr${i}.txt \
  -o ${wrkdir}/ibdmix_out/gtex_ibdmix_chr${i}.txt \
  -w &
done

# convert genotype files to bed files
for i in {1..22}
do
  (cat ${wrkdir}/generate_gt_out/altai_gtex_combined_gt_chr${i}.txt |\
  sed -e1,1d |\
  awk 'BEGIN { FS="\t"; OFS="\t" } { $2=$2 "\t" $2 } 1' \
  > ${wrkdir}/tag_snps/altai_gtex_combined_gt_chr${i}.bed ;
  bgzip ${wrkdir}/tag_snps/altai_gtex_combined_gt_chr${i}.bed ;
  tabix -p bed ${wrkdir}/tag_snps/altai_gtex_combined_gt_chr${i}.bed.gz) &
done

mkdir -p ${wrkdir}/tag_snps/tag_snp_query/
mkdir -p ${wrkdir}/tag_snps/tag_snp_gt/

# compute tag SNP stats with R script





