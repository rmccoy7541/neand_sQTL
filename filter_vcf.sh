# This script assumes that the genotype file has already been downloaded. It filters non-biallelic sites and indexes it.

## Step 4 - Genotype & Covariates Preparation
################################################
# insert scripts dir below
scripts=$(echo /work-zfs/rmccoy22/aseyedi2/neanderthal-sqtl/src/primary/)
# insert the full path for GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz
vcf=$(echo /scratch/groups/rmccoy22/Ne_sQTL/files/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz)
# WGS directory
wgs=$(echo /scratch/groups/rmccoy22/Ne_sQTL/files/)

# filter the non-biallelic sites from genotype file using bcftools; see script for details
sbatch --wait --export=vcf=$vcf ${scripts}/sh/bcftools_filter.sh
# index our friend with tabix
echo "Indexing our friend..."
interact -p shared -t 360:00 -c 6
ml htslib; tabix -p vcf $wgs/GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz
exit
