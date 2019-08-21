# module load python/3.7
# python -m venv ./venv
# source ./venv/bin/activate
# pip install snakemake
# snakemake --version

configfile: "config.yaml"

rule all:
    input: 
        expand("Ne-sQTL_perind.counts.gz.qqnorm_chr{i}.qtltools",i=range(1,22))

rule filter_vcf:
    input:
        vcf=expand("{vcf}", vcf=config["vcf"]),
        ncbiFiles=expand("{ncbiFiles}", ncbiFiles=config["ncbiFiles"])
    output:
        expand("{ncbiFiles}/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTExWGSGenotypeMatrixBiallelicOnly.HQ.vcf.gz", ncbiFiles=config["ncbiFiles"])
    threads: 23 # in addition to the 1 thread, so 24 total
    shell:
        "bcftools view -m2 -M2 -v snps --threads {threads} -O z -o {output} {input.vcf}"

rule index_vcf:
    input:
        expand("{ncbiFiles}/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTExWGSGenotypeMatrixBiallelicOnly.HQ.vcf.gz", ncbiFiles=config["ncbiFiles"])
    output:
        expand("{ncbiFiles}/phg000830.v1.GTEx_WGS.genotype-calls-vcf.c1/GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz.tbi",ncbiFiles=config["ncbiFiles"]),
        touch(expand("{ncbiFiles}/.index_vcf.chkpnt",ncbiFiles=config["ncbiFiles"]))
    shell:
        "tabix -p vcf {input}"


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
    message:
        "Intron clustering..."
    params:
        LC=config["leafcutter"]
    shell:
        "mkdir intronclustering/;"
        "python {params.LC}/clustering/leafcutter_cluster.py \
        -j juncfiles.txt \
        -r intronclustering/ \
        -m 50 \
        -o Ne-sQTL \
        -l 500000"

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
        "python {params.LC}/scripts/prepare_phenotype_table.py Ne-sQTL_perind.counts.gz -p 10"

rule QTLtools_filter:
    input:
         file=expand("Ne-sQTL_perind.counts.gz.qqnorm_chr{i}",i=range(1,22))
    output:
        expand("{input.file}.qtltools")
    message:
        "Making phenotype files QTLtools compatible..."
    shell:
        "cat {input.file} | awk '{ $4=$4\" . +\"; print $0 }' | tr " " \"\t\" | bgzip -c > {input.file}.qtltools"

rule index_phen:
    input:
        expand("Ne-sQTL_perind.counts.gz.qqnorm_chr{i}.qtltools",i=range(1,22))
    output:
        expand("Ne-sQTL_perind.counts.gz.qqnorm_chr{i}.qtltools.tbi",i=range(1,22))
    message:
        "Indexing phenotype files..."
    shell:
        "tabix -p bed {input}"

rule sra_tissue_xtract:
    input:
        "metadata/SraRunTable.txt",
        "metadata/GTExTissueKey.csv"
    output:
        "tissue_table.txt"
    message:
        "Extracting tested tissues..."
    script:
        "src/sqtl_mapping/primary/R/05_sraTissueExtract.R"

#######
## Starting from the end
#######
