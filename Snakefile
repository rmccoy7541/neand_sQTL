configfile: "config.yaml"
#localrules:

rule filter_vcf:
    input:
        vcf=config["vcf"],
        ncbiFiles=config["ncbiFiles"]
    output:
        "{ncbiFiles}GTExWGSGenotypeMatrixBiallelicOnly.HQ.vcf.gz"
    shell:
        "sbatch --wait --export=vcf={vcf},outdir=$PWD src/sqtl_mapping/primary/sh/00a_bcftools_filter.sh"

rule index_vcf:
    input:
        "GTExWGSGenotypeMatrixBiallelicOnly.HQ.vcf.gz"
    output:
        "GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz.tbi"
    shell:
        "sbatch --export=outdir=$PWD src/sqtl_mapping/primary/sh/00b_index_vcf.sh"

rule junc_cluster:
    input:
        LC=config["leafcutter"]
    output:
        ".junc_cluster.chkpnt"
    shell:
        "sbatch --wait src/sqtl_mapping/sh/01_junc_cluster.sh"

rule intron_clustering:
    input:
        LC=config["leafcutter"]
    output:
        ".intron_clustering.chkpnt"
    shell:
        "sbatch --wait src/sqtl_mapping/sh/02_intronclustering.sh;"
        "cd intronclustering/"

rule prepare_phen_table:
    input:
        LC=config["leafcutter"]
    output:
        ".prepare_phen_table.chkpnt"
    shell:
        "sbatch --wait /sh/03_prepare_phen_table.sh {LC}"






