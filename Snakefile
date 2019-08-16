# module load python/3.7
# python -m venv ./venv
# source ./venv/bin/activate
# pip install snakemake
# snakemake --version

configfile: "config.yaml"
#localrules:

vcf=config["vcf"],
ncbiFiles=config["ncbiFiles"]
LC=config["leafcutter"]

rule all:
     input: ".prepare_phen_table.chkpnt"

rule filter_vcf:
    output:
        "{ncbiFiles}/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTExWGSGenotypeMatrixBiallelicOnly.HQ.vcf.gz"
    threads:
        24
    shell:
        "bcftools view -m2 -M2 -v snps --threads 23 -O z -o {output} {vcf}"

rule index_vcf:
     input:
        "{ncbiFiles}/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTExWGSGenotypeMatrixBiallelicOnly.HQ.vcf.gz"
     output:
        "{ncbiFiles}/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz.tbi",
        touch("{ncbiFiles}/.index_vcf.chkpnt")
     shell:
        "tabix -p vcf {input}"
#
# rule junc_cluster:
#     input:
#         expand("{ncbiFiles}/.index_vcf.chkpnt", ncbiFiles=config["ncbiFiles"])
#     output:
#         ".junc_cluster.chkpnt"
#     shell:
#         "sbatch --wait src/sqtl_mapping/sh/01_junc_cluster.sh;"
#         "touch .junc_cluster.chkpnt"
#
# rule intron_clustering:
#     input:
#         ".junc_cluster.chkpnt"
#     output:
#         ".intron_clustering.chkpnt"
#     shell:
#         expand("sbatch --wait src/sqtl_mapping/sh/02_intronclustering.sh {LC};"
#         "touch .intron_clustering.chkpnt;"
#         "cd intronclustering/", LC=config["leafcutter"])
#
# rule prepare_phen_table:
#     input:
#         LC,
#         ".intron_clustering.chkpnt"
#     output:
#         ".prepare_phen_table.chkpnt"
#     shell:
#         expand("sbatch --wait src/sqtl_mapping/sh/03_prepare_phen_table.sh {LC};"
#         "touch .prepare_phen_table.chkpnt",LC=config["leafcutter"])
#
#
#
#
#

