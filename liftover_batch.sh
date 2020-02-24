#!/bin/bash
#SBATCH --job-name=LiftoverVcf
#SBATCH --time=12:0:0
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=1
#SBATCH --array=1-22

ml java

java -jar ~/work/progs/gatk-4.0.12.0/gatk-package-4.0.12.0-local.jar LiftoverVcf \
  -I /scratch/groups/rmccoy22/aseyedi2/NL_sQTL_iso/IBDmix/archaicGenomes/vcfChrom/AltaiNea.hg19_1000g.$SLURM_ARRAY_TASK_ID.LowQualRemoved.vcf \
  -O /scratch/groups/rmccoy22/aseyedi2/NL_sQTL_iso/IBDmix/archaicGenomes/AltaiNea.hg38_1000g.$SLURM_ARRAY_TASK_ID.LowQualRemoved.liftover.vcf \
  -C /scratch/groups/rmccoy22/aseyedi2/NL_sQTL_iso/IBDmix/archaicGenomes/b37ToHg38.over.chain \
  -R /scratch/groups/rmccoy22/aseyedi2/NL_sQTL_iso/IBDmix/archaicGenomes/hg38.fa \
  --REJECT /scratch/groups/rmccoy22/aseyedi2/NL_sQTL_iso/IBDmix/archaicGenomes/rejectedRecords_chr$SLURM_ARRAY_TASK_ID.vcf \
  --ALLOW_MISSING_FIELDS_IN_HEADER \
  --TMP_DIR /scratch/groups/rmccoy22/aseyedi2/NL_sQTL_iso/IBDmix/archaicGenomes/tempdir_$SLURM_ARRAY_TASK_ID \
  &> liftOver_wTempDir_chr$SLURM_ARRAY_TASK_ID.log
