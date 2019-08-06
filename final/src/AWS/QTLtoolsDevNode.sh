#!/bin/bash

ml htslib
ml R
ml gcc

# this project's scripts dir
scripts=$(echo /home-1/aseyedi2@jhu.edu/work/aseyedi2/neand_sQTL/src/primary/)
# data dir
data=$(echo /home-1/aseyedi2@jhu.edu/work/aseyedi2/neand_sQTL/data/)
# ncbi/files/
ncbiFiles=$(echo /scratch/groups/rmccoy22/Ne_sQTL/files/)
# IF YOU ALREADY HAVE NON-BIALLELIC INDEXED VCF
VCF=$(echo /scratch/groups/rmccoy22/Ne_sQTL/files/GTExWGSGenotypeMatrixBiallelicOnly.vcf.gz)
# input directory with junc files here
junc=$(echo '/scratch/groups/rmccoy22/Ne_sQTL/sra/sqtl_junc')
# leafcutter directory here
leafCutter=$(echo /scratch/groups/rmccoy22/aseyedi2/leafcutter)


# so strange... I had to use `bind 'set disable-completion on'` to make it work
{
for line in $(cat GTExCovKey.csv); do
	full=$(echo $line | awk -F',' '{print $1}')
	abb=$(echo $line | awk -F',' '{print $2}')
	if grep "$abb" tissuesused.txt; then
		for i in {1..100}; do
			/scratch/groups/rmccoy22/progs/QTLtools/QTLtools_1.1_Ubuntu14.04_x86_64 cis \
			  --vcf  $VCF \
			  --bed $abb/$abb.pheno.bed.gz \
			  --cov  $abb/$full.v7.covariates_output.txt \
			  --nominal 1  \
			  --chunk $i 100 \
			  --out "$abb/${abb}_nominals_chunk_$i.txt"
		done

		for i in {1..100}; do
			/home-1/aseyedi2@jhu.edu/work/progs/QTLtools/QTLtools_1.1_Ubuntu14.04_x86_64 cis \
			  --vcf  $VCF \
			  --bed $abb/$abb.pheno.bed.gz \
			  --cov  $abb/$full.v7.covariates_output.txt \
			  --permute 1000  \
			  --chunk $i 100 \
			  --out "$abb/${abb}_permutations_chunk_$i.txt" \
			  --normal
	  	done

		Rscript ~/work/progs/QTLtools/script/runFDR_cis.R $abb/$abb.permutations_full.txt.gz 0.05 $abb/$abb.permutations_full_FDR

		for i in {1..100}; do
			/home-1/aseyedi2@jhu.edu/work/progs/QTLtools/QTLtools_1.1_Ubuntu14.04_x86_64 cis \
			  --vcf  $VCF \
			  --bed $abb/$abb.pheno.bed.gz \
			  --cov  $abb/$full.v7.covariates_output.txt \
			  --mapping $abb/$abb.permutations_full_FDR.thresholds.txt \
			  --chunk $i 100 \
			  --out "$abb/${abb}_conditionals_chunk_$i.txt"
		done

		# catting noms	
		for i in {1..100}; do
			cat $abb/${abb}_nominals_chunk_${i}.txt | gzip -c >> $abb/$abb.nominals.all.chunks.txt.gz
		done

		# getting chunks
		ls $abb/${abb}_* | sort -V > $abb/${abb}_chunks_list.txt
		#######

		for i in {1..100}; do
			cat ${abb}/${abb}_nominals_chunk_${i}_out.txt | gzip -c >> ${abb}/$abb.nominals.all.chunks.NE_only.txt.gz
		done

		mkdir $abb/nominals; mv $abb/*_nominals_* $abb/nominals/

		#######

		for i in {1..100}; do
			cat $abb/${abb}_permutations_chunk_$i.txt | gzip -c >> $abb/${abb}.permutations_full.txt.gz
		done

		mkdir ${abb}/permutations; mv ${abb}/*_permutations_* ${abb}/permutations/

		#######

		for i in {1..100}; do
			cat $abb/${abb}_conditionals_chunk_$i.txt | gzip -c >> $abb/${abb}.conditional_full.txt.gz
		done

		# Consider passing WHLBLD_chunks as a command line var
		for line in $(cat $abb/${abb}_chunks_list.txt); do
			Rscript ${scripts}/R/NomPassExtract.R $abb/${line} "${data}/../analysis/SPRIME/sprime_calls.txt"
		done

		mkdir /home-1/aseyedi2@jhu.edu/work/aseyedi2/sQTL/$abb

		Rscript ${scripts}/R/QQPlot-Viz.R /home-1/aseyedi2@jhu.edu/work/aseyedi2/sQTL/$abb $abb/$abb.nominals.all.chunks.NE_only.txt.gz $abb/$abb.permutations_full.txt.gz ${data}/../analysis/SPRIME/sprime_calls.txt

		else
			echo "$abb not found..."
	fi
done
} > QTLTools.log.04.22.2019.AWStest.txt
