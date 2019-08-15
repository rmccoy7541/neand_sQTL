configfile: "config.yaml"
#localrules:
vcf=config["vcf"],
ncbiFiles=config["ncbiFiles"]
LC=config["leafcutter"]

rule all:
    input: ".prepare_phen_table.chkpnt"

rule filter_vcf:
    input:
        vcf,
        ncbiFiles
    output:
        "{ncbiFiles}/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTExWGSGenotypeMatrixBiallelicOnly.HQ.vcf.gz"
    shell:
        "sbatch --wait --export=vcf={vcf},outdir=$PWD src/sqtl_mapping/primary/sh/00a_bcftools_filter.sh"

rule index_vcf:
    input:
        "{ncbiFiles}/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTExWGSGenotypeMatrixBiallelicOnly.HQ.vcf.gz"
    output:
        "{ncbiFiles}/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz.tbi",
        "{ncbiFiles}/.index_vcf.chkpnt"
    shell:
        "sbatch --export=outdir=$PWD src/sqtl_mapping/primary/sh/00b_index_vcf.sh;"
        "{ncbiFiles}/touch .index_vcf.chkpnt"

rule junc_cluster:
    input:
        "{ncbiFiles}/.index_vcf.chkpnt"
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
        "sbatch --wait src/sqtl_mapping/sh/02_intronclustering.sh {LC};"
        "touch .intron_clustering.chkpnt;"
        "cd intronclustering/"

rule prepare_phen_table:
    input:
        LC
    output:
        ".prepare_phen_table.chkpnt"
    shell:
        "sbatch --wait src/sqtl_mapping/sh/03_prepare_phen_table.sh {LC};"
        "touch .prepare_phen_table.chkpnt"






