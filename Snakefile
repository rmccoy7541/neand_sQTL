# module load python/3.7
# python -m venv ./venv
# source ./venv/bin/activate
# pip install snakemake
# snakemake --version

configfile: "config.yaml"

rule all:
    input: 
        "leafcutterphenotypes.txt"

rule filter_vcf:
    input:
        vcf=expand("{vcf}", vcf=config["vcf"]),
        ncbiFiles=expand("{ncbiFiles}", ncbiFiles=config["ncbiFiles"])
    output:
        expand("{ncbiFiles}/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTExWGSGenotypeMatrixBiallelicOnly.HQ.vcf.gz", ncbiFiles=config["ncbiFiles"])
    threads: 23 # in addition to the 1 thread, so 24 total
    shell:
         #ml bcftools
        "bcftools view -m2 -M2 -v snps --threads {threads} -O z -o {output} {input.vcf}"

rule index_vcf:
    input:
        expand("{ncbiFiles}/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTExWGSGenotypeMatrixBiallelicOnly.HQ.vcf.gz", ncbiFiles=config["ncbiFiles"])
    output:
        expand("{ncbiFiles}/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz.tbi",ncbiFiles=config["ncbiFiles"]),
        touch(expand("{ncbiFiles}/.index_vcf.chkpnt",ncbiFiles=config["ncbiFiles"]))
    shell:
        #ml htslib
        "tabix -p vcf $outdir/GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz"


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
    message:
        "Preparing phenotype table..."
    params:
        LC=config["leafcutter"]
    shell:
        "src/sqtl_mapping/sh/03_prepare_phen_table.sh {params.LC}"

rule write_LC_phen:
    input:
        ".prepare_phen_table.chkpnt"
    output:
        "leafcutterphenotypes.txt"
    shell:
        "ls *qqnorm* > {output}"

