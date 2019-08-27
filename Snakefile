# module load python/3.7
# python -m venv ./venv
# source ./venv/bin/activate
# pip install snakemake
# snakemake --version

def read_tissues_output(wildcard):
    with open(wildcard) as f:
        samples = [sample for sample in f.read().split('\n') if len(sample) > 0]  # we dont want empty lines
        return samples

# snakemake --dag -n -F | dot -Tsvg > dag.svg

configfile: "config.yaml"

rule all:
    input:
        ".sort_zip_ind_pheno.chkpnt"

# not sure how this rule below works with the inputs and outputs being so vague but ok
rule QTLtools_filter:
    input:
        phen="{i}.gz",
        chk=".prepare_phen_table.chkpnt"
    output:
        "{i}.qtltools"
    message:
        "Making phenotype files QTLtools compatible..."
    shell:
        "cat {input.phen} | awk '{{ $4=$4\" . +\"; print $0 }}' | tr " " \"\t\" | bgzip -c > {output}"

rule sra_tissue_xtract:
    input:
        "metadata/SraRunTable.txt",
        "metadata/GTExTissueKey.csv",
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
        "tissue_table.txt",
    output:
        "tissuenames.txt"
    shell:
        "cat {input} | cut -f3 | awk '{{ if(NR>1)print }}' |  awk '!seen[$0]++' > {output}"

rule make_tis_dirs:
    input:
        ".sra_name_change.chkpnt",

    output:
        touch(".make_tis_dirs.chkpnt")
    message:
        "Making directories for each type of tissue, saving tissue types in a file, and "
        "moving each outputted file into its respective tissue folder..."
    shell: # maybe figure out regex for 1_*.txt? I feel like that's inefficient
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

rule cat_phen_files:
    input:
        tis=read_tissues_output('tissuesused.txt'),
        chk=".move_tis.chkpnt"
    output:
        touch(".cat_phen_files.chkpnt")
    shell: # I have no idea if this will work lol - also figure out how to parallelize this
        "head -1 {input.tis}/1_{input.tis}.txt > $line/$line.phen_fastqtl.bed;"
        "for file in $(ls {input.tis}/*_*.txt); do"
            "cat $file | sed -e1,1d >> {input.tis}/{input.tis}.phen_fastqtl.bed;"
        "done"

rule sort_zip_ind_pheno:
    input:
        tis=read_tissues_output('tissuesused.txt'),
        chk2=".cat_phen_files.chkpnt"
    output:
        touch(".sort_zip_ind_pheno.chkpnt")
    shell:
        "bedtools sort -header -i {input.tis}/{input.tis}.phen_fastqtl.bed > \
        {input.tis}/{input.tis}.pheno.bed;"
        "bgzip -f {input.tis}/{input.tis}.pheno.bed;"
        "tabix -p bed {input.tis}/{input.tis}.pheno.bed.gz;"
        "mkdir {input.tis}/sepfiles;"
        "mv {input.tis}/*_{input.tis} {input.tis}/sepfiles/"

rule DL_cov:
    input:
        cov="https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_sQTL_covariates.tar.gz"
    output:
        "GTEx_Analysis_v8_sQTL_covariates.tar.gz"
    shell:
        "wget {input.cov}"

rule move_cov:
    input:
        cov=read_tissues_output('metadata/GTExCovKey.csv'),
        unpack="GTEx_Analysis_v8_sQTL_covariates.tar.gz",
        tis="tissuesused.txt"
    output:
        touch(".move_cove.chkpnt")
    shell:
        "tar -xvf {input.unpack};"
        "full=$(echo {input.cov} | awk -F',' '{{ print $1 }}');"
        "abb=$(echo {input.cov} | awk -F',' '{{ print $2 }}');"
        "if grep \"$abb\" {input.tis}; then"
        "cp GTEx_Analysis_v8_sQTL_covariates/$full.v7.covariates.txt $abb;"
        "done"

rule QTLtools_Loop:
    input:
        vcf=expand("{vcf}", vcf=config["vcf"])
    output:
        touch(".QTLtools_loop.chkpnt")
    shell: # include bash script but parallelize alone two different dimensions
        "sbatch --wait \
        --export=scripts=$scripts,data=$data,vcf=$vcf,sprime=$sprime \
        -a 2-$numTissues \
        ${scripts}/sh/08_QTLTools-Loop.sh"



# for line in $(cat tissuesused.txt)
# do
#    mkdir $line/sepfiles
#    mv $line/*_${line}.txt $line/sepfiles/
# done
#
# # download genotype covariates
# wget https://storage.googleapis.com/gtex_analysis_v7/single_tissue_eqtl_data/GTEx_Analysis_v7_eQTL_covariates.tar.gz
#
# tar -xvf GTEx_Analysis_v7_eQTL_covariates.tar.gz
#
# cp $data/GTExCovKey.csv $PWD
#
# # Moves covariates to corresponding directory
# for line in $(cat GTExCovKey.csv)
# do
#    full=$(echo $line | awk -F',' '{print $1}')
#    abb=$(echo $line | awk -F',' '{print $2}')
#    if grep "$abb" tissuesused.txt; then
#       cp GTEx_Analysis_v7_eQTL_covariates/$full.v7.covariates.txt $abb
#       Rscript ${scripts}/R/07_mergePCs.R Ne-sQTL_perind.counts.gz.PCs $abb/$full.v7.covariates.txt tissuetable/tissue_table.txt
#       mv $full.v7.covariates_output.txt $abb
#    fi
# done
#
# ## Step 4 - Mapping sQTLs using QTLtools
# ################################################
# numTissues=$(wc -l GTExCovKey.csv)
#
# # Will take at least 3 weeks lol
# sbatch --wait \
#   --export=scripts=$scripts,data=$data,vcf=$vcf,sprime=$sprime \
#   -a 2-$numTissues \
#   ${scripts}/sh/08_QTLTools-Loop.sh
#
# mkdir -p $backupdir
# mkdir $backupdir/all_noms
# mkdir $backupdir/sqtl_nrich
# mv $PWD/*/*permutation* $backupdir
# mv $PWD/*/*pheno* $backupdir
# mv tissuenames.txt $backupdir
#
# cd $backupdir
# cp $sprime .
#
# # concatenates permutation pass results
# for TISSUE in ADPSBQ ADPVSC ADRNLG ARTACRN ARTAORT ARTTBL BREAST BRNACC BRNAMY BRNCDT BRNCHA BRNCHB BRNCTXA BRNCTXB BRNHPP BRNHPT BRNNCC BRNPTM BRNSNG BRNSPC CLNSGM CLNTRN ESPGEJ ESPMCS ESPMSL FIBRLBLS HRTAA HRTLV LCL LIVER LUNG MSCLSK NERVET OVARY PNCREAS PRSTTE PTTARY SKINNS SKINS SLVRYG SNTTRM SPLEEN STMACH TESTIS THYROID UTERUS VAGINA WHLBLD
# do
#   mkdir ${TISSUE}
#   mv ${TISSUE}*permutations* ${TISSUE}
#   for i in {1..100}
#   do
#     echo "Catting ${TISSUE}_permutations_chunk_${i}.txt..."
#     cat ${TISSUE}/${TISSUE}_permutations_chunk_${i}.txt >> ${TISSUE}_permutations.txt
#   done
# done
#
# # getting tissue names
# for i in $(ls *_permutations.txt | sort -V); do echo $i | cut -d'_' -f 1; done > tissuenames.txt
#
# wget ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Assembly/GRCh37-HG19_Broad_variant/Homo_sapiens_assembly19.fasta
#
# # get allele frequences from VCF
# sbatch --export=VCF=$VCF $scripts/sh/13_VariantToTable.sh
#
# # concat the GTEx AF VCF chunks
# cat GTExWGS.AF.chr1.txt > GTExWGS.AF.all.txt
# for i in {2..22}; do tail -n +2 GTExWGS.AF.chr${i}.txt >> GTExWGS.AF.all.txt ; done
# rm GTExWGS.AF.chr*.txt

#######
## Starting from the end
#######
