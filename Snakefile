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
    input:
        expand("{vcf}", vcf=config["vcf"]),
        expand("{ncbiFiles}", ncbiFiles=config["ncbiFiles"])
    output:
        expand("{ncbiFiles}/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTExWGSGenotypeMatrixBiallelicOnly.HQ.vcf.gz", ncbiFiles=config["ncbiFiles"])
    shell:
        expand("sbatch --wait --export=vcf={vcf},outdir=$PWD src/sqtl_mapping/primary/sh/00a_bcftools_filter.sh", vcf=config["vcf"])

rule index_vcf:
    input:
        expand("{ncbiFiles}/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTExWGSGenotypeMatrixBiallelicOnly.HQ.vcf.gz", ncbiFiles=config["ncbiFiles"])
    output:
        expand("{ncbiFiles}/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz.tbi",ncbiFiles=config["ncbiFiles"]),
        expand("{ncbiFiles}/.index_vcf.chkpnt",ncbiFiles=config["ncbiFiles"])
    shell:
        expand("sbatch --export=outdir=$PWD src/sqtl_mapping/primary/sh/00b_index_vcf.sh;"
        "touch {ncbiFiles}/.index_vcf.chkpnt",ncbiFiles=config["ncbiFiles"])

rule junc_cluster:
    input:
        expand("{ncbiFiles}/.index_vcf.chkpnt", ncbiFiles=config["ncbiFiles"])
    output:
        ".junc_cluster.chkpnt"
    shell:
        "sbatch --wait src/sqtl_mapping/sh/01_junc_cluster.sh;"
        "touch .junc_cluster.chkpnt"

rule intron_clustering:
    input:
        ".junc_cluster.chkpnt"
    output:
        ".intron_clustering.chkpnt"
    shell:
        expand("sbatch --wait src/sqtl_mapping/sh/02_intronclustering.sh {LC};"
        "touch .intron_clustering.chkpnt;"
        "cd intronclustering/", LC=config["leafcutter"])

rule prepare_phen_table:
    input:
        LC,
        ".intron_clustering.chkpnt"
    output:
        ".prepare_phen_table.chkpnt"
    shell:
        expand("sbatch --wait src/sqtl_mapping/sh/03_prepare_phen_table.sh {LC};"
        "touch .prepare_phen_table.chkpnt",LC=config["leafcutter"])






