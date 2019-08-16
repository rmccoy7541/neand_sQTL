configfile: "config.yaml"

rule all:
    input: 
        ".prepare_phen_table.chkpnt"

rule filter_vcf:
    input:
        expand("{vcf}", vcf=config["vcf"]),
        expand("{ncbiFiles}", ncbiFiles=config["ncbiFiles"])
    output:
        expand("{ncbiFiles}/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTExWGSGenotypeMatrixBiallelicOnly.HQ.vcf.gz", ncbiFiles=config["ncbiFiles"])
    shell:
        "src/sqtl_mapping/primary/sh/00a_bcftools_filter.sh"


rule index_vcf:
    input:
        expand("{ncbiFiles}/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTExWGSGenotypeMatrixBiallelicOnly.HQ.vcf.gz", ncbiFiles=config["ncbiFiles"])
    output:
        expand("{ncbiFiles}/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz.tbi",ncbiFiles=config["ncbiFiles"]),
        touch(expand("{ncbiFiles}/.index_vcf.chkpnt",ncbiFiles=config["ncbiFiles"]))
    shell:
        "src/sqtl_mapping/primary/sh/00b_index_vcf.sh"


rule junc_cluster:
    input:
        expand("{ncbiFiles}/.index_vcf.chkpnt", ncbiFiles=config["ncbiFiles"])
    output:
        touch(".junc_cluster.chkpnt")
    shell:
        "src/sqtl_mapping/sh/01_junc_cluster.sh"


rule intron_clustering:
    input:
        ".junc_cluster.chkpnt"
    output:
        touch(".intron_clustering.chkpnt")
    params:
        LC=config["leafcutter"]
    shell:
        "src/sqtl_mapping/sh/02_intronclustering.sh {params.LC}"


rule prepare_phen_table:
    input:
        config["leafcutter"],
        ".intron_clustering.chkpnt"
    output:
        touch(".prepare_phen_table.chkpnt")
    params:
        LC=config["leafcutter"]
    shell:
        "src/sqtl_mapping/sh/03_prepare_phen_table.sh {params.LC}"