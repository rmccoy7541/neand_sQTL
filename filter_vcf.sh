# This script assumes that the genotype file has already been downloaded. It filters non-biallelic sites and indexes it.

## Step 4 - Genotype & Covariates Preparation
################################################
scripts=$(echo /home-1/aseyedi2@jhu.edu/work/aseyedi2/neand_sQTL/src/primary/)

# filter the non-biallelic sites from genotype file using bcftools; see script for details
sbatch --wait ${scripts}/sh/bcftools_filter.sh /scratch/groups/rmccoy22/Ne_sQTL/files
# index our friend with tabix
echo "Indexing our friend..."
interact -p shared -t 360:00 -c 6
ml htslib; tabix -p vcf /scratch/groups/rmccoy22/Ne_sQTL/files/GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz
exit
