# module load python/3.7
# python -m venv ./venv
# source ./venv/bin/activate
# pip install snakemake
# snakemake --version

# snakemake --dag -n | dot -Tsvg > dag.svg

configfile: "config.yaml"

include: "sprime_run_Snakemake"

rule all:
    input:
        "sprime_calls.txt",
        "GTEx_Analysis_v8_sQTL/",
        "GTEx_Analysis_v8_sQTL_phenotype_matrices/"
#     input:
#         "GTEx_Analysis_v8_sQTL/",
#         "GTEx_Analysis_v8_sQTL_phenotype_matrices/",
#         expand("kg_vcf/1kg_yri_chr{q}.vcf.gz", q=range(1,23)),
#         "kg_vcf/1kg_yri_chrX.vcf.gz",
#         expand("{kg_dir}/ALL.chr{q}.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz", kg_dir=config["kg_dir"], q=range(1,23)),
#         expand("gtex_vcf/gtex_chr{v}.snps.recode.vcf.gz", v=range(1,23)),
#         "gtex_vcf/gtex_chrX.snps.recode.vcf.gz.tbi"


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

rule over_chain:
    output:
        "metadata/hg38ToHg19.over.chain"
    shell:
        "wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz metadata/;"
        "unzip {output}"

rule sprime_R:
    input:
        results=rules.sprime_run.output,
        arch_vcf=config["arch_vcf"],
        over_chain=rules.over_chain.output
    output:
        "sprime_calls.txt"
    script:
        "src/sprime/sprime_neand.R"
