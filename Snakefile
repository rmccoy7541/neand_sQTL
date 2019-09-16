# module load python/3.7
# python -m venv ./venv
# source ./venv/bin/activate
# pip install snakemake
# snakemake --version

# snakemake --dag -n | dot -Tsvg > dag.svg

configfile: "config.yaml"

include: "sprime_run_Snakemake"

TISSUES = ["Adipose_Subcutaneous",
"Adipose_Visceral_Omentum",
"Adrenal_Gland",
"Artery_Aorta",
"Artery_Coronary",
"Artery_Tibial",
"Brain_Amygdala",
"Brain_Anterior_cingulate_cortex_BA24",
"Brain_Caudate_basal_ganglia",
"Brain_Cerebellar_Hemisphere",
"Brain_Cerebellum",
"Brain_Cortex",
"Brain_Frontal_Cortex_BA9",
"Brain_Hippocampus",
"Brain_Hypothalamus",
"Brain_Nucleus_accumbens_basal_ganglia",
"Brain_Putamen_basal_ganglia",
"Brain_Spinal_cord_cervical_c-1",
"Brain_Substantia_nigra",
"Breast_Mammary_Tissue",
"Cells_Cultured_fibroblasts",
"Cells_EBV-transformed_lymphocytes",
"Colon_Sigmoid",
"Colon_Transverse",
"Esophagus_Gastroesophageal_Junction",
"Esophagus_Mucosa",
"Esophagus_Muscularis",
"Heart_Atrial_Appendage",
"Heart_Left_Ventricle",
"Kidney_Cortex",
"Liver",
"Lung",
"Minor_Salivary_Gland",
"Muscle_Skeletal",
"Nerve_Tibial",
"Ovary",
"Pancreas",
"Pituitary",
"Prostate",
"Skin_Not_Sun_Exposed_Suprapubic",
"Skin_Sun_Exposed_Lower_leg",
"Small_Intestine_Terminal_Ileum",
"Spleen",
"Stomach",
"Testis",
"Thyroid",
"Uterus",
"Vagina",
"Whole_Blood"]

rule all:
    input:
        "sprime_calls.txt",
        "GTEx_Analysis_v8_sQTL/",
        "GTEx_Analysis_v8_sQTL_phenotype_matrices/"
        # expand("{tissue}_permutation_table_NE.txt", tissue=TISSUES)
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
        "GTEx_Analysis_v8_sQTL_phenotype_matrices/",
        touch(".decomp.chkpnt")
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

# rule neand_sQTL:
#     input:
#         chk=".decomp.chkpnt",
#         perm=expand("GTEx_Analysis_v8_sQTL/{tissue}.v8.sqtl_signifpairs.txt.gz", tissue=TISSUES),
#         sprime="sprime_calls.txt"
#     output:
#         "{tissue}_permutation_table_NE.txt"
#     script:
#         "src/analysis/NE_sQTL.R"

