# module load python/3.7
# python -m venv ./venv
# source ./venv/bin/activate
# pip install snakemake
# snakemake --version

# snakemake --dag -n | dot -Tsvg > dag.svg

configfile: "config.yaml"

include: "sprime_run_Snakemake"

TISSUES = ["Adipose_Subcutaneous", # 763
"Adipose_Visceral_Omentum", # 564
"Adrenal_Gland", # 275
"Artery_Aorta", # 450
"Artery_Coronary", # 253
"Artery_Tibial", # 770
"Brain_Amygdala", # 177
"Brain_Anterior_cingulate_cortex_BA24", # 213
"Brain_Caudate_basal_ganglia", # 291
"Brain_Cerebellar_Hemisphere", # 263
"Brain_Cerebellum", # 298
"Brain_Cortex", # 325
"Brain_Frontal_Cortex_BA9", # 425
"Brain_Hippocampus", # 243
"Brain_Hypothalamus", # 236
"Brain_Nucleus_accumbens_basal_ganglia", # 277
"Brain_Putamen_basal_ganglia", # 232
"Brain_Spinal_cord_cervical_c-1", # 182
"Brain_Substantia_nigra", # 164
"Breast_Mammary_Tissue", # 480
"Cells_Cultured_fibroblasts", # 527
"Cells_EBV-transformed_lymphocytes", # 192
"Colon_Sigmoid", # 389
"Colon_Transverse", # 432
"Esophagus_Gastroesophageal_Junction", # 401
"Esophagus_Mucosa", # 622
"Esophagus_Muscularis", # 559
"Heart_Atrial_Appendage", # 452
"Heart_Left_Ventricle", # 689
"Kidney_Cortex", # 100
"Liver", # 251
"Lung", # 867
"Minor_Salivary_Gland", # 181
"Muscle_Skeletal", # 1132
"Nerve_Tibial", # 722
"Ovary", # 195
"Pancreas", # 360
"Pituitary", # 301
"Prostate", # 262
"Skin_Not_Sun_Exposed_Suprapubic", # 638
"Skin_Sun_Exposed_Lower_leg", # 849
"Small_Intestine_Terminal_Ileum", # 193
"Spleen", # 260
"Stomach", # 381
"Testis", # 406
"Thyroid", # 812
"Uterus", # 166
"Vagina", # 173
"Whole_Blood"] # 3288

rule all:
    input:
        "sQTLs_per_tissue.png",
        "TopGenes_PermPass_All.csv"

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
        "metadata/sprime_calls.txt"
    script:
        "src/sprime/sprime_neand.R"

rule dl_gtf:
    output:
        "gencode.v26.GRCh38.genes.gtf"
    params:
        url="https://storage.googleapis.com/gtex_analysis_v8/reference/gencode.v26.GRCh38.genes.gtf"
    shell:
        "wget {params.url}"

rule neand_sQTL:
    input:
        chk=".decomp.chkpnt",
        perm=expand("GTEx_Analysis_v8_sQTL/{tissue}.v8.sqtl_signifpairs.txt.gz", tissue=TISSUES),
        sprime="metadata/sprime_calls.txt",
        gtf="gencode.v26.GRCh38.genes.gtf"
    output:
        "{tissue}_permutation_table_NE.txt"
    script:
        "src/analysis/NE_sQTL.R"

rule count_sQTL:
    input:
        expand("{tissue}_permutation_table_NE.txt", tissue=TISSUES)
    output:
        "sQTLs_per_tissue.png",
        "TopGenes_PermPass_All.csv"
    script:
        "src/analysis/count_sqtl.R"

# TODO
# rule manhattan:
#     input:
#         expand("{tissue}_permutation_table_NE.txt", tissue=TISSUES),
#         "metadata/sprime_calls.txt"
#     output:

rule dl_intronCounts:
    output:
        "GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz"
    params:
        introns="https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz"
    shell:
        "wget {params.introns}"

rule process_introns:
    input:
        dl_intronCounts.output
    output:
        "GTEx_v8_junctions_nohead.gct.gz"
    shell:
        "zcat {input} | grep -v \"#1.2\" | grep -v '357746'$'\t''17382 > GTEx_v8_junctions_nohead.gct;"
        "bgzip GTEx_v8_junctions_nohead.gct"

rule splitIntronCounts:
    input:
        introns=process_introns.output,
        tistab="metadata/tissue_key.csv"
    output:
        "{tissue}_intronCounts.txt"
    script:
        "src/analysis/preprocess_intronCounts.R"

rule find_NL_introns:
    input:
        introns=expand("{tissue}_intronCounts.txt", tissue=TISSUES),
        # TODO vcf= Find out how steph preprocessed the VCF
        # 