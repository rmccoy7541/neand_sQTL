#!/bin/bash

cd /work-zfs/rmccoy22/rmccoy22/sqtl/enrichment/tissue_vcf

for TISSUE in ADPSBQ ADPVSC ADRNLG ARTACRN ARTAORT ARTTBL BREAST BRNACC BRNAMY BRNCDT BRNCHA BRNCHB BRNCTXA BRNCTXB BRNHPP BRNHPT BRNNCC BRNPTM BRNSNG BRNSPC CLNSGM CLNTRN ESPGEJ ESPMCS ESPMSL FIBRLBLS HRTAA HRTLV LCL LIVER LUNG MSCLSK NERVET OVARY PNCREAS PRSTTE PTTARY SKINNS SKINS SLVRYG SNTTRM SPLEEN STMACH TESTIS THYROID UTERUS VAGINA WHLBLD
do
echo ${TISSUE}
cat /work-zfs/rmccoy22/aseyedi2/sqtl_permutation_backup/${TISSUE}_permutations.txt | cut -f8 -d ' ' | sort -k1,1n -k2,2n -t "_" > ../snp_lists/${TISSUE}.list
/work-zfs/rmccoy22/progs/gatk-4.1.2.0/gatk \
  SelectVariants \
  --variant /scratch/groups/rmccoy22/Ne_sQTL/files/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz \
  --keep-ids ../snp_lists/${TISSUE}.list \
  --sample-name ../sample_lists/${TISSUE}.args \
  --allow-nonoverlapping-command-line-samples \
  --output ${TISSUE}.vcf &
done
