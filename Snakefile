# module load python/3.7
# python -m venv ./venv
# source ./venv/bin/activate
# pip install snakemake
# snakemake --version

configfile: "config.yaml"

rule all:
    input:
        expand("gtex_vcf/gtex_chr{i}.vcf",i=range(1,22)),
        "gtex_vcf/gtex_chrX.vcf",
        expand("kg_vcf/1kg_yri_chr{q}.vcf.gz", q=range(1,22)),
        "kg_vcf/1kg_yri_chrX.vcf.gz",
        # expand("{kg_dir}/ALL.chr{p}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz", kg_dir=config["kg_dir"], p=range(1,22))

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

### Consider making a sub workflow

# rule mkdir_vcf:
#     output:
#         directory("gtex_vcf/"),
#         directory("kg_vcf/")
#     shell:
#         "mkdir -p {output}"

rule vcf_split1_23:
    input:
        vcf=config["vcf"]
    output:
        "gtex_vcf/gtex_chr{i}.vcf"
    threads:
        23
    shell:
        "tabix -h {input.vcf} chr{wildcards.i} > {output}"

rule YRI_select:
    input:
        "metadata/yri.txt"
    output:
        "kg_vcf/1kg_yri_chr{q}.vcf.gz",
        # "{kg_dir}/ALL.chr{p}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz"
    shell:
        "bcftools view --force-samples "
        "-S {input} "
        "-V indels "
        "-O z "
        "-o kg_vcf/1kg_yri_chr{wildcards.q}.vcf.gz "
        # "{wildcards.kg_dir}/ALL.chr{wildcards.p}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz"
