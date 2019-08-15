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
chkpntDir=config["chkpntDir"]



rule all:
    input: ".prepare_phen_table.chkpnt"

rule filter_vcf:
    input:
        vcf,
        ncbiFiles
    output:
        "{ncbiFiles}/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTExWGSGenotypeMatrixBiallelicOnly.HQ.vcf.gz"
    shell:
        "mkdir -p {chkpntDir};"
        "sbatch --wait --export=vcf={vcf},outdir=$PWD src/sqtl_mapping/primary/sh/00a_bcftools_filter.sh"

rule index_vcf:
    input:
        "{ncbiFiles}/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTExWGSGenotypeMatrixBiallelicOnly.HQ.vcf.gz"
    output:
        "{ncbiFiles}/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz.tbi",
        "{chkpntDir}/.index_vcf.chkpnt"
    shell:
        "sbatch --export=outdir=$PWD src/sqtl_mapping/primary/sh/00b_index_vcf.sh;"
        "{chkpntDir}/touch .index_vcf.chkpnt"

rule junc_cluster:
    input:
        "{chkpntDir}/.index_vcf.chkpnt"
    output:
        "{chkpntDir}/.junc_cluster.chkpnt"
    shell:
        "sbatch --wait src/sqtl_mapping/sh/01_junc_cluster.sh;"
        "touch {ncbiFiles}/.junc_cluster.chkpnt"

rule intron_clustering:
    input:
        "{chkpntDir}/.junc_cluster.chkpnt"
    output:
        "{chkpntDir}/.intron_clustering.chkpnt"
    shell:
        "sbatch --wait src/sqtl_mapping/sh/02_intronclustering.sh {LC};"
        "touch {chkpntDir}/.intron_clustering.chkpnt;"
        "cd intronclustering/"

rule prepare_phen_table:
    input:
        LC,
        "{chkpntDir}/.intron_clustering.chkpnt"
    output:
        "{chkpntDir}/.prepare_phen_table.chkpnt"
    shell:
        "sbatch --wait src/sqtl_mapping/sh/03_prepare_phen_table.sh {LC};"
        "touch {chkpntDir}/.prepare_phen_table.chkpnt"






