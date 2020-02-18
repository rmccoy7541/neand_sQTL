ml bcftools

for i in {1..22}; do
  wget "http://cdna.eva.mpg.de/neandertal/altai/AltaiNeandertal/VCF/AltaiNea.hg19_1000g.${i}.mod.vcf.gz"
done

for i in {1..22}; do
  gunzip AltaiNea.hg19_1000g.${i}.mod.vcf.gz
done

vcf-concat -f concatList.txt > AltaiNea.hg19_100g.concat.vcf

# now liftover the vcf
wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/liftOver/hg19ToHg38.over.chain.gz

sed '/LowQual/d' ./AltaiNea.hg38_1000g.concat.vcf > AltaiNea.hg19_1000g.concat.LowQualRemoved.vcf

java -jar ~/work/progs/gatk-4.0.12.0/gatk-package-4.0.12.0-local.jar LiftoverVcf \
  -I AltaiNea.hg19_1000g.concat.LowQualRemoved.vcf \
  -O AltaiNea.hg38_1000g.concat.LowQualRemoved.vcf \
  -C b37ToHg38.over.chain \
  -R hg38.fa \
  --REJECT rejectedRecords.vcf \
  --ALLOW_MISSING_FIELDS_IN_HEADER \
  --TMP_DIR /scratch/groups/rmccoy22/aseyedi2/NL_sQTL_iso/IBDmix/archaicGenomes/tempdir \
  &> liftOver_wTempDir.log

bcftools view -m 2 -M 2 -v snps --threads 23 -O v -o GTExWGSGenotypeMatrixBiallelicOnly_v8.vcf.gz /scratch/groups/rmccoy22/aseyedi2/sQTLv8/data/vcf/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz

# From IBDMix
generate_gt -a AltaiNea.hg38_1000g.concat.LowQualRemoved.vcf -m GTExWGSGenotypeMatrixBiallelicOnly_v8.vcf -o IBDMix_GT.txt