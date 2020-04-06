# module load python/3.7
# python -m venv ./venv
# source ./venv/bin/activate
# pip install snakemake
# snakemake --version

# snakemake --dag -n | dot -Tsvg > dag.svg

# Force
# snakemake --dag -n -Ft | dot -Tsvg > dag.svg

configfile: "config.yaml"

include: "sprime_run_Snakemake"
include: "preprocessVCF"

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
        "TopGenes_PermPass_All.csv",
        "metadata/sprime_calls.txt",
        "results/SeparateTissues.png",
        "results/AllTissuesNLIso.png",
        "results/loosenedRestrictionsGenes.txt"

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

# Result of Snakemake Pipeline
rule sprime_R:
    input:
        results=expand("{sprime_dir}/output/results.chr{z}.score", sprime_dir=config["sprime_dir"], z=range(1,23)),
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
        perm=expand("{sQTLs}/{tissue}.v8.sqtl_signifpairs.txt.gz", sQTLs=config["sQTLs"], tissue=TISSUES),
        sprime="metadata/sprime_calls.txt",
        gtf="gencode.v26.GRCh38.genes.gtf"
    output:
        "results/perms/{tissue}_permutation_table_NE.txt"
    script:
        "src/analysis/NE_sQTL.R"

rule count_sQTL:
    input:
        expand("results/perms/{tissue}_permutation_table_NE.txt", tissue=TISSUES)
    output:
        "results/sQTLs_per_tissue.png",
        "results/TopGenes_PermPass_All.csv"
    script:
        "src/analysis/count_sqtl.R"

# 1. Change header of GTEx VCF - shadow rule? place the output somewhere else?
rule get_vcf_header:
    input:
        vcf=config["vcf"]
    output:
        "new_gtex_header.hdr"
    shell:
        "bcftools view -h {input} > {output};"

rule change_vcf_header:
    input:
        header=rules.get_vcf_header.output,
        vcf=rules.get_vcf_header.input
    output:
        "GTEx_v8_reheader.vcf.gz"
    shell:
        "bcftools reheader -h {input.header} -o {output} {input.vcf}"

rule tabix_new_vcf:
    input:
        rules.change_vcf_header.output
    output:
        "GTEx_v8_reheader.vcf.gz.tbi"
    shell:
        "tabix -p vcf {input}"

# 2. Subset VCF to contain only variants that were called as introgressed by sprime
rule extract_SPrime:
    input:
        "metadata/sprime_calls.txt"
    output:
        "sprime_variants.list"
    shell:
        "cat {input} | grep -v \"CHROM\" | cut -f 4 > {output}"

rule dl_intronCounts:
    output:
        "GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz"
    params:
        introns="https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz"
    shell:
        "wget {params.introns}"

rule process_introns:
    input:
        "GTEx_Analysis_2017-06-05_v8_STARv2.5.3a_junctions.gct.gz"
    output:
        "GTEx_v8_junctions_nohead.gct.gz"
    shell:
        "zcat {input} | grep -v \"#1.2\" | grep -v '357746'$'\t''17382 > GTEx_v8_junctions_nohead.gct;"
        "bgzip GTEx_v8_junctions_nohead.gct"

rule splitIntronCounts:
    input:
        introns="GTEx_v8_junctions_nohead.gct.gz",
        tistab="metadata/tissue_key.csv"
    output:
        "results/IC/{tissue}_intronCounts.txt"
    script:
        "src/analysis/preprocess_intronCounts.R"

rule find_NL_introns:
    input:
        introns=expand("results/IC/{tissue}_intronCounts.txt", tiss ue=TISSUES),
        perm=expand("{sQTLs}/{tissue}.v8.sqtl_signifpairs.txt.gz", tissue=TISSUES, sQTLs=config["sQTLs"]),
        vcf_merge="vcf_for_merge.txt.gz",
    output:
        "results/finalIsos/{tissue}_NL_isos.txt"
    script:
        "src/analysis/merge_tables.R"

rule countCounts:
    input:
         finalISO=expand("results/finalIsos/{tissue}_NL_isos.txt", tissue=TISSUES)
    output:
        "results/countCounts/{tissues}_countCounts.txt"
    script:
        "src/analysis/countCounts.R"

rule plotCountCounts:
    input:
        countCountsDir="results/countCounts/"
    output:
        "results/SeparateTissues.png",
        "results/AllTissuesNLIso.png"
    script:
        "src/analysis/plot_countCounts.R"

rule getRelScripts:
    input:
        wd="results/countCounts/"
    output:
        "results/loosenedRestrictions.txt"
    script:
        "src/analysis/getRelevantTranscripts.R"

rule findGenes:
    input:
        relScripts=rules.getRelScripts.output,
        gtf="gencode.v26.GRCh38.genes.gtf"
    output:
        "results/loosenedRestrictionsGenes.txt"
    script:
        "src/analysis/findGenes.R"

