#!/usr/bin/env bash

# note: currently written to be executed on development node
# could modify as an array job for SLURM

ml bcftools
ml vcftools
ml htslib
ml java

cd /scratch/users/rmccoy22@jhu.edu/gtex_sprime

# get a list of Yoruban samples from 1kg project (negligible Neanderthal introgression) 
cat /scratch/groups/rmccoy22/resources/1kg-hs37d5/integrated_call_samples_v3.20130502.ALL.panel | grep 'YRI' | cut -f1 > yri.txt

# split the GTEx VCF by chromosome
for i in {1..22}
do
  tabix GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz ${i} > gtex_chr${i}.vcf &
done

wait

# select only SNPs, YRI samples from 1kg VCFs
for i in {1..22}
do
vcftools \
  --gzvcf /scratch/groups/rmccoy22/resources/1kg-hs37d5/ALL.chr${i}.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz \
  --keep yri.txt \
  --remove-indels \
  --recode \
  --out 1kg_yri_chr${i} &
done

wait

# select only SNPs from GTEx VCFs
for i in {1..22}
do
vcftools \
  --vcf gtex_chr${i}.vcf \
  --remove-indels \
  --recode \
  --out gtex_chr${i} &
done

wait

# bgzip and tabix index
for i in *.recode.vcf
do
  bgzip ${i} && tabix -p vcf ${i}.gz &
done

wait

# merge 1kg_YRI and GTEx
for i in {1..22}
do
  bcftools filter \
    --missing-to-ref \
    -O z \
    -o merged_unfiltered_chr${i}.vcf.gz \
    1kg_yri_chr${i}.recode.vcf.gz \
    gtex_chr${i}.recode.vcf.gz &
done

wait

# filter (keep PASS only)
for i in {1..22}
do
  bcftools view \
    --apply-filters .,PASS \
    -O z \
    -o merged_filtered_chr${i}.vcf.gz \
    merged_unfiltered_chr${i}.vcf.gz &
done

wait

# concatenate genetic maps
for i in {1..22}
do 
  cat genetic_map/plink.chr${i}.GRCh37.map >> genetic_map/plink.all_autosomes.GRCh37.map
done

# get a list of merged, filtered VCFs to concatenate
for i in {1..22}
do
  echo "/scratch/users/rmccoy22@jhu.edu/archaic_splicing/gtex_sprime/merged_filtered_chr"${i}".vcf.gz" >> processed_vcf_list.txt
done

# concatenate VCFs
bcftools concat \
  --file-list processed_vcf_list.txt \
  --naive \
  --threads 24 \
  -O z \
  -o merged_filtered_all_autosomes.vcf.gz

# index
tabix -p vcf merged_filtered_all_autosomes.vcf.gz

# exclude SNPs with missing genotypes and exclude SNPs with AC=0 (no ALT genotypes)
# apply phase 3 pilot mask from 1000 Genomes, indicating that the site was accessible to genotyping
bcftools filter \
  --include 'AN=1520 && AC > 0' \
  --regions-file 20140520.pilot_mask.autosomes.bed \
  --threads 24 \
  -O z \
  -o merged_filtered_all_autosomes_an1520.vcf.gz \
  merged_filtered_all_autosomes.vcf.gz

# run sprime
for i in {1..22}
do   
  java -Xmx12g -jar ~/work/progs/sprime.jar \
    gt=merged_filtered_all_autosomes_an1520.vcf.gz \
    chrom=${i} \
    outgroup=yri.txt \
    map=genetic_map/plink.all_autosomes.GRCh37.map \
    out=output/results.chr${i} & 
done

wait

echo "SPrime complete."