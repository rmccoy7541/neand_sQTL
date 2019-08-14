#!/bin/bash
#SBATCH --partition=shared
#SBATCH --job-name=QTLTools-Int
#SBATCH --time=24:0:0
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=2
# enter array info command line

line=`sed "${SLURM_ARRAY_TASK_ID}q;d" GTExCovKey.csv`

full=$(echo $line | awk -F',' '{print $1}')
abb=$(echo $line | awk -F',' '{print $2}')

# catting noms
for i in {1..100}; do
 cat $abb/${abb}_nominals_chunk_${i}.txt | gzip -c >> $abb/$abb.nominals.all.chunks.txt.gz
done

# # getting chunks
# # ls $abb/${abb}_* | sort -V >> $abb/${abb}_chunks_list.txt
# #######

# for i in {1..100}; do
# cat ${abb}/${abb}_nominals_chunk_${i}_out.txt | gzip -c >> ${abb}/$abb.nominals.all.chunks.NE_only.txt.gz
# done

# # mkdir $abb/nominals; mv $abb/*_nominals_* $abb/nominals/

#######

for i in {1..100}; do
  cat $abb/${abb}_permutations_chunk_$i.txt | gzip -c >> $abb/${abb}.permutations_full.txt.gz
done

mkdir ${abb}/permutations; mv ${abb}/*_permutations_* ${abb}/permutations/

#######

# for i in {1..100}; do
  # cat $abb/${abb}_conditionals_chunk_$i.txt | gzip -c >> $abb/${abb}.conditional_full.txt.gz
# done

# mkdir ${abb}/conditionals; mv ${abb}/*_conditionals_* ${abb}conditionals/

mkdir /home-1/aseyedi2@jhu.edu/work/aseyedi2/sQTL/$abb
