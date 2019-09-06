# module load python/3.7
# python -m venv ./venv
# source ./venv/bin/activate
# pip install snakemake
# snakemake --version

configfile: "config.yaml"

subworkflow sprime_prep:
    workdir:
           "src/sprime"
    snakefile:
            "src/sprime/Snakefile"
    configfile:
            "src/sprime/sprime_config.yaml"

rule all:
    input:
        "GTEx_Analysis_v8_sQTL/",
        "GTEx_Analysis_v8_sQTL_phenotype_matrices/",



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
#         directory("kg_vcf/"),
#         touch(".mkdir_vcf.chkpnt")
#     shell:
#         "mkdir -p {output}"


