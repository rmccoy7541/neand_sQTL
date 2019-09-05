# module load python/3.7
# python -m venv ./venv
# source ./venv/bin/activate
# pip install snakemake
# snakemake --version

configfile: "config.yaml"

# gcloud auth application-default login

rule all:
    input:
        expand("gtex_vcf/gtex_chr{i}.vcf",i=range(1,22)),
        "gtex_vcf/gtex_chrX.vcf"

rule dl_files:
    params:
        sqtl="https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_sQTL.tar",
        phen="https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_sQTL_phenotype_matrices.tar"
    output:
        "GTEx_Analysis_v8_sQTL.tar",
        "GTEx_Analysis_v8_sQTL_phenotype_matrices.tar"
    shell:
        "wget {params.sqtl}; wget {params.phen}"

rule decomp:
    input:
        sqtl="GTEx_Analysis_v8_sQTL.tar",
        phen="GTEx_Analysis_v8_sQTL_phenotype_matrices.tar"
    output:
        "GTEx_Analysis_v8_sQTL/",
        "GTEx_Analysis_v8_sQTL_phenotype_matrices/"
    shell:
        "tar -xvf {input.sqtl}; tar -xvf {input.phen};"
        "rm {input.sqtl}; rm {input.phen}"

rule mkdir_vcf:
    output:
        "gtex_vcf/",
        "kg_vcf/",
        touch(".mkdir_vcf.chkpnt")
    shell:
        "mkdir -p {output}"

rule vcf_split1_23:
    input:
        vcf=config["vcf"]
    output:
        "gtex_vcf/gtex_chr{i}.vcf"
    threads:
        23
    shell:
        "tabix -h {input.vcf} chr{wildcards.i} > {output}"

