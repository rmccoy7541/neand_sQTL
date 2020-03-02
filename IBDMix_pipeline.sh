ml bcftools
ml plink
ml java

for i in {1..22}; do
  wget "http://cdna.eva.mpg.de/neandertal/altai/AltaiNeandertal/VCF/AltaiNea.hg19_1000g.${i}.mod.vcf.gz"
done

for i in {1..22}; do
  gunzip AltaiNea.hg19_1000g.${i}.mod.vcf.gz
done

#don't need to concat
#vcf-concat -f concatList.txt > AltaiNea.hg19_1000g.concat.vcf

# now liftover the vcf
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz

for i in {1..22}; do
  sed '/LowQual/d' ./AltaiNea.hg19_1000g.${i}.mod.vcf > AltaiNea.hg19_1000g.${i}.LowQualRemoved.vcf
done

for i in {1..22}; do
  java -jar ~/work/progs/gatk-4.0.12.0/gatk-package-4.0.12.0-local.jar LiftoverVcf \
    -I AltaiNea.hg19_1000g.${i}.LowQualRemoved.vcf \
    -O AltaiNea.hg38_1000g.${i}.LowQualRemoved.vcf \
    -C ../b37ToHg38.over.chain \
    -R ../hg38.fa \
    --REJECT rejectedRecords_chr${i}.vcf \
    --ALLOW_MISSING_FIELDS_IN_HEADER \
    --TMP_DIR /scratch/groups/rmccoy22/aseyedi2/NL_sQTL_iso/IBDmix/archaicGenomes/tempdir \
    &> liftOver_wTempDir_chr${i}.log
done

bcftools view -m 2 -M 2 -v snps --threads 23 -O v -o GTExWGSGenotypeMatrixBiallelicOnly_v8.vcf.gz /scratch/groups/rmccoy22/aseyedi2/sQTLv8/data/vcf/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz

# From IBDMix
# Try taking the first MB of chromosome 1 and see if that works with generate_gt
# maybe strip it down to essential parts of VCF - maybe remove metadata
# make a bed file of the first MB and do bedtools intersect on both those files to just extract that bit
# Consider doing it on individual chromosomes separately
plink --vcf GTExWGSGenotypeMatrixBiallelicOnly_v8.vcf --recode vcf --out gtex_plink_recode

for i in {1..22}; do
  plink --vcf AltaiNea.hg38_1000g.${i}.LowQualRemoved.vcf --recode vcf --out altai_plink_recode_chr${i} --allow-extra-chr
done

for i in {1..22}; do
  ../../generate_gt -a altai_plink_recode_chr${i}.vcf  -m ../gtex_plink_recode.vcf -o ../IBDMixOut/IBDMix_GT_chr${i}.vcf
done
