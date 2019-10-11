#!/bin/bash
#SBATCH --job-name=sprime_run
#SBATCH --time=24:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=4
#SBATCH --array=1-22

### 02-sprime-run.sh ###

i=${SLURM_ARRAY_TASK_ID}

ml bcftools
ml vcftools
ml htslib
ml java

# gtex-sprime directory
wrkdir=/work-zfs/rmccoy22/aseyedi2/sprime_v7_test
# 1kg project direcotry (1kg-hs37d5/)
kg_dir=/work-zfs/rmccoy22/resources/1kg-hs37d5/
# original gtex vcf (prefilter) GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz
gtex_vcf=/scratch/groups/rmccoy22/aseyedi2/sQTLv8/ncbi/dbGaP-20712/files/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz

cd ${wrkdir}

mkdir -p output

java -Xmx12g -jar ~/work/progs/sprime.jar \
  gt=merged_filtered_all_autosomes.vcf.gz \
  chrom=${i} \
  outgroup=yri.txt \
  map=genetic_map/plink.all_autosomes.GRCh37.map \
  out=output/results.chr${i}

touch sprime_v7_run_complete