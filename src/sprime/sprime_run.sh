#!/bin/bash
#SBATCH --job-name=sprime_run
#SBATCH --time=24:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=4
#SBATCH --array=1-22%22

### 02-sprime-run.sh ###

if [ ${SLURM_ARRAY_TASK_ID} -eq 23 ]
then
    i="X"
else
    i=${SLURM_ARRAY_TASK_ID}
fi

ml bcftools
ml vcftools
ml htslib
ml java

# gtex-sprime directory
wrkdir=/scratch/groups/rmccoy22/aseyedi2/sQTLv8/sprime
# 1kg project direcotry (1kg-hs37d5/)
kg_dir=/work-zfs/rmccoy22/resources/1kg-GRCh38
# original gtex vcf (prefilter) GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_635Ind_PASS_AB02_GQ20_HETX_MISS15_PLINKQC.vcf.gz
gtex_vcf=/scratch/groups/rmccoy22/aseyedi2/sQTLv8/data/vcf/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.vcf.gz

cd ${wrkdir}

mkdir -p output

java -Xmx12g -jar ~/work/progs/sprime.jar \
  gt=merged_filtered_all_autosomes.vcf.gz \
  chrom=chr${i} \
  outgroup=yri.txt \
  map=genetic_map/plink.all_autosomes.GRCh38.map \
  out=output/results.chr${i}


# java -Xmx12g -jar ~/work/progs/sprime.jar \
#   gt=/scratch/groups/rmccoy22/aseyedi2/sQTLv8/sprime/merged_filtered_all_autosomes.vcf.gz \
#   outgroup=/scratch/groups/rmccoy22/aseyedi2/sQTLv8/sprime/yri.txt \
#   map=/scratch/groups/rmccoy22/aseyedi2/sQTLv8/sprime/genetic_map/plink.all_autosomes.GRCh38.map \
#   out=/scratch/groups/rmccoy22/aseyedi2/sQTLv8/sprime_test/results.test.ADAMSTL3.minscore \
#   chrom=chr15:83537709-84537709 \
#   minscore=50000