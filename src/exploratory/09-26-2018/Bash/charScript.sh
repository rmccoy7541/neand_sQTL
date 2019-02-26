#!/bin/sh
#SBATCH ^^partition=shared
#SBATCH --time=00:360:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mail-type=end
#SBATCH --mail-user=aseyedi2@jhu.edu

zcat ~/work/1kg/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '#' | cut -f 2-3 > ~/data/aseyedi2/1kg_data/s01/rsID_chr1.txt
zcat ~/work/1kg/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '#' | cut -f 2-3 > ~/data/aseyedi2/1kg_data/s02/rsID_chr2.txt
zcat ~/work/1kg/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '#' | cut -f 2-3 > ~/data/aseyedi2/1kg_data/s03/rsID_chr3.txt
zcat ~/work/1kg/ALL.chr4.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '#' | cut -f 2-3 > ~/data/aseyedi2/1kg_data/s04/rsID_chr4.txt
zcat ~/work/1kg/ALL.chr5.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '#' | cut -f 2-3 > ~/data/aseyedi2/1kg_data/s05/rsID_chr5.txt
zcat ~/work/1kg/ALL.chr6.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '#' | cut -f 2-3 > ~/data/aseyedi2/1kg_data/s06/rsID_chr6.txt
zcat ~/work/1kg/ALL.chr7.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '#' | cut -f 2-3 > ~/data/aseyedi2/1kg_data/s07/rsID_chr7.txt
zcat ~/work/1kg/ALL.chr8.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '#' | cut -f 2-3 > ~/data/aseyedi2/1kg_data/s08/rsID_chr8.txt
zcat ~/work/1kg/ALL.chr9.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '#' | cut -f 2-3 > ~/data/aseyedi2/1kg_data/s09/rsID_chr9.txt
zcat ~/work/1kg/ALL.chr10.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '#' | cut -f 2-3 > ~/data/aseyedi2/1kg_data/s10/rsID_chr10.txt
zcat ~/work/1kg/ALL.chr11.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '#' | cut -f 2-3 > ~/data/aseyedi2/1kg_data/s11/rsID_chr11.txt
zcat ~/work/1kg/ALL.chr12.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '#' | cut -f 2-3 > ~/data/aseyedi2/1kg_data/s12/rsID_chr12.txt
zcat ~/work/1kg/ALL.chr13.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '#' | cut -f 2-3 > ~/data/aseyedi2/1kg_data/s13/rsID_chr13.txt
zcat ~/work/1kg/ALL.chr14.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '#' | cut -f 2-3 > ~/data/aseyedi2/1kg_data/s14/rsID_chr14.txt
zcat ~/work/1kg/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '#' | cut -f 2-3 > ~/data/aseyedi2/1kg_data/s15/rsID_chr15.txt
zcat ~/work/1kg/ALL.chr16.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '#' | cut -f 2-3 > ~/data/aseyedi2/1kg_data/s16/rsID_chr16.txt
zcat ~/work/1kg/ALL.chr17.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '#' | cut -f 2-3 > ~/data/aseyedi2/1kg_data/s17/rsID_chr17.txt
zcat ~/work/1kg/ALL.chr18.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '#' | cut -f 2-3 > ~/data/aseyedi2/1kg_data/s18/rsID_chr18.txt
zcat ~/work/1kg/ALL.chr19.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '#' | cut -f 2-3 > ~/data/aseyedi2/1kg_data/s19/rsID_chr19.txt
zcat ~/work/1kg/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '#' | cut -f 2-3 > ~/data/aseyedi2/1kg_data/s20/rsID_chr20.txt
zcat ~/work/1kg/ALL.chr21.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '#' | cut -f 2-3 > ~/data/aseyedi2/1kg_data/s21/rsID_chr21.txt
zcat ~/work/1kg/ALL.chr22.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz | grep -v '#' | cut -f 2-3 > ~/data/aseyedi2/1kg_data/s22/rsID_chr22.txt
