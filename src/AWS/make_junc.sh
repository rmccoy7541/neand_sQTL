#!/bin/bash
#$ -t 1-11688
#$ -tc 500

export PATH="/root/progs/leafcutter/scripts:$PATH"
export PATH="/root/progs/samtools-1.9:$PATH"

echo "Task id is ${SGE_TASK_ID}."

SRR_ID="`sed -n ${SGE_TASK_ID}p reference/srr_list.txt`"

fusera mount --token ~/keys/prj_20712.ngc --accession "${SRR_ID}" ~/gtex_v7 &

yes "" | sleep 5s

cram_input=`basename gtex_v7/${SRR_ID}/*.cram .cram`

bam2junc.sh gtex_v7/${SRR_ID}/${cram_input}.cram ${cram_input}.junc reference/Homo_sapiens_assembly19.fasta

awk '{if ($1 !~ "^GL" && $1 !~ "^NC") print}' ${cram_input}.junc > ${cram_input}_filtered.junc

yes | scp -i "keys/mykey.rsa" "${cram_input}_filtered.junc" root@ec2-54-227-207-203.compute-1.amazonaws.com:~/complete/

fusera unmount gtex_v7