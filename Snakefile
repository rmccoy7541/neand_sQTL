# module load python/3.7
# python -m venv ./venv
# source ./venv/bin/activate
# pip install snakemake
# snakemake --version

def fileAsList(file):
    with open(file) as f:
        for line in f:
            lis = []
            spl = line.split()
            lis.append(spl[0])
        return lis


configfile: "config.yaml"

rule all:
    input: 
        ".move_tis.chkpnt"

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
        touch(".prepare_phen_table.chkpnt"),
        expand("Ne-sQTL_perind.counts.gz.qqnorm_chr{i}.gz",i=range(1,22))
    message:
        "Preparing phenotype table..."
    params:
        LC=config["leafcutter"]
    shell:
        "python {params.LC}/scripts/prepare_phenotype_table.py Ne-sQTL_perind.counts.gz -p 10"

rule QTLtools_filter:
    input:
        "{i}.gz",
        ".prepare_phen_table.chkpnt"
    output:
        "{i}.qtltools"
    message:
        "Making phenotype files QTLtools compatible..."
    shell:
        "cat {input} | awk '{{ $4=$4\" . +\"; print $0 }}' | tr " " \"\t\" | bgzip -c > {output}"

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

rule sra_name_change_sort:
    input:
        expand("Ne-sQTL_perind.counts.gz.qqnorm_chr{i}.qtltools",i=range(1,22)),
        "tissue_table.txt"
    output:
        touch(".sra_name_change.chkpnt")
    script:
        "src/sqtl_mapping/primary/R/06_sraNameChangeSort.R"

rule get_tis_names:
    input:
        "tissue_table.txt"
    output:
        "tissuenames.txt"
    shell:
        "cat {input} | cut -f3 | awk '{{if(NR>1)print}}' |  awk '!seen[$0]++' > {output}"

rule make_tis_dirs:
    input:
        ".sra_name_change.chkpnt"
    output:
        touch(".make_tis_dirs.chkpnt")
    message:
        "Making directories for each type of tissue, saving tissue types in a file, and "
        "moving each outputted file into its respective tissue folder..."
    shell:
        "for i in 1_*.txt; do echo $i | cut -d'_' -f 2| cut -d'.' -f 1 | xargs mkdir; done;"

checkpoint get_tissue:
    input:
        ".sra_name_change.chkpnt"
    output:
        "tissuesused.txt"
    shell:
        "for i in 1_*.txt; do echo $i | cut -d'_' -f 2| cut -d'.' -f 1 >> {output}; done;"

rule move_tis:
    input:
        ".make_tis_dirs.chkpnt"
    output:
        touch(".move_tis.chkpnt")
    shell:
        "for i in *_*.txt; do echo $i | awk -F'[_.]' '{{print $2}}' | xargs -I '{{}}' mv $i '{{}}' ; done"

# def read_tissues_output():
#     with open('tissuesused.txt') as f:
#         samples = [sample for sample in f.read().split('\n') if len(sample) > 0]  # we dont want empty lines
#         return samples
#
# rule sort_zip_ind_pheno:
#     input:
#         read_tissues_output(),
#         ".move_tis.chkpnt"
#     output:
#         touch(".sort_zip_ind_pheno.chkpnt")
#     shell:
#         "bedtools sort -header -i {input.tis}/{input.tis}.phen_fastqtl.bed > \
#         {input.tis}/{input.tis}.pheno.bed;"
#         "bgzip -f {input.tis}/{input.tis}.pheno.bed;"
#         "tabix -p bed {input.tis}/{input.tis}.pheno.bed.gz"




#######
## Starting from the end
#######
