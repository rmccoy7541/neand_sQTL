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
        "GTEx_Analysis_v8_sQTL_phenotype_matrices/",
        # rules.index_merged.output,
        # rules.cat_genetic_maps.output,
        # expand("{sprime_dir}/output/results.chr{z}.score", sprime_dir=config["sprime_dir"], z=range(1,23))

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
    params:
        url="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz"
    shell:
        "wget {params.url} metadata/;"
        "gunzip {output}.gz"

rule sprime_R:
    input:
        results=expand("{sprime_dir}/output/results.chr{z}.score", sprime_dir=config["sprime_dir"], z=range(1,23)),
        arch_vcf=config["arch_vcf"],
        over_chain=rules.over_chain.output
    output:
        "sprime_calls.txt"
    script:
        "src/sprime/sprime_neand.R"
