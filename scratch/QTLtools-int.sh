### this performs all of the intermediary steps of QTLtools


# catting noms
for i in {1..100}; do
cat $abb/${abb}_nominals_chunk_${i}.txt | gzip -c >> $abb/$abb.nominals.all.chunks.txt.gz
done

# getting chunks
ls $abb/${abb}_* | sort -V >> $abb/${abb}_chunks_list.txt
#######

rm $abb/sprime_calls.txt

for i in {1..100}; do
cat ${abb}/${abb}_nominals_chunk_${i}_out.txt | gzip -c >> ${abb}/$abb.nominals.all.chunks.NE_only.txt.gz
done

mkdir $abb/nominals; mv $abb/*_nominals_* $abb/nominals/

#######

for i in {1..100}; do
  cat $abb/${abb}_permutations_chunk_$i.txt | gzip -c >> $abb/${abb}.permutations_full.txt.gz
done

mkdir ${abb}/permutations; mv ${abb}/*_permutations_* ${abb}/permutations/