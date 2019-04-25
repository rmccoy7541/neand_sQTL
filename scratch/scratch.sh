sbatch --wait --export=VCF=$VCF,pheno=$(echo BRNCHA.pheno.bed.gz),tissue=$(echo BRNCHA),covariates=$(echo Brain_Cerebellum.v7.covariates_output.txt) ${scripts}/sh/NomPass.sh

sbatch --export=listPath=$PWD,tissue=$(echo BRNCHA),scripts=$scripts ${scripts}sh/NomPassExtractCall.sh

Rscript ${scripts}/R/NomPassExtract.R TESTIS.nominals.all.chunks.txt "sprime_calls.txt"

`split -d -a 6 -l 65898466 --additional-suffix '.split' TESTIS.nominals.all.chunks.txt ''`

sbatch --export=listPath=$PWD,tissue=$(echo TESTIS),scripts=$scripts ../NomPassExtractCall.sh

for i in {1..100}; do
	cat BRNCHA_nominals_chunk_${i}_out.txt | gzip -c >> BRNCHA.nominals.all.chunks.NE_only.txt.gz
done

for i in $(cat ~/failedpermpass); do find -iname slurm-$i.out | xargs head -n 50 | grep "Chunk\|phenotype data in\|Read sample list from\|Reading sample list from" > out; chunk=$(cat out | grep "Chunk" | awk -F'=' '{print $2}' | cut -d/ -f1 | cut -d[ -f2); pheno=$(cat out | grep "Reading sample list" | awk -F"[" '{print $2}' | awk -F"]" '{print $1}'| grep "pheno"); cov=$(cat out|grep "covariates"| cut -d[ -f2 | cut -d] -f1); sbatch -a $SLURM_ARRAY_TASK_ID --export=VCF=$VCF,pheno=$(echo $abb/$abb.pheno.bed.gz),tissue=$(echo $abb/$abb),covariates=$(echo $abb/$full.v7.covariates_output.txt),abb=$abb,full=$full ${scripts}/sh/PermPass.sh; done



for i in $(cat ~/failedpermpass); do find -iname slurm-$i.out | xargs head -n 50 | grep "Chunk\|phenotype data in\|Read sample list from\|Reading sample list from" > out; chunk=$(cat out | grep "Chunk" | awk -F'=' '{print $2}' | cut -d/ -f1 | cut -d[ -f2); pheno=$(cat out | grep "Reading sample list" | awk -F"[" '{print $2}' | awk -F"]" '{print $1}'| grep "pheno"); cov=$(cat out|grep "covariates"| cut -d[ -f2 | cut -d] -f1); tissue=$(echo $pheno | cut -d. -f1); sbatch -a $chunk --export=VCF=$VCF,pheno=$(echo $pheno),tissue=$(echo $tissue),covariates=$(echo $cov) ${scripts}/sh/PermPass.sh; done