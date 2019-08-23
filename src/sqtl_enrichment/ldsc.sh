#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=ldsc
#SBATCH --nodes=1
#SBATCH --time=12:0:0
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --array=1-22

ml htslib
ml plink
ml python/2.7-anaconda

CHROM=${SLURM_ARRAY_TASK_ID}

# use tabix to extract chromosome from VCF
tabix -h \
/scratch/groups/rmccoy22/Ne_sQTL/files/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_652Ind_GATK_HaplotypeCaller.vcf.gz \
${CHROM} \
> /work-zfs/rmccoy22/rmccoy22/sqtl/enrichment/ldsc/chr${CHROM}.vcf

plink \
--cm-map /work-zfs/rmccoy22/resources/reference/genetic_maps/chr${CHROM}.b37.gmap ${CHROM} \
--make-bed \
--out chr${CHROM}_w_cms \
--vcf chr${CHROM}.vcf

conda activate ldsc

python /work-zfs/rmccoy22/progs/ldsc/ldsc.py \
--bfile /work-zfs/rmccoy22/rmccoy22/sqtl/enrichment/ldsc/chr${CHROM}_w_cms \
--l2 \
--ld-wind-cm 1 \
--out /work-zfs/rmccoy22/rmccoy22/sqtl/enrichment/ldsc/chr${CHROM}_w_cms

