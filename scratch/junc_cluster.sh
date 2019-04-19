#!/bin/bash
for i in $(ls *junc); do
	new=$(echo $i | awk -F'.' '{print $1}')
	mv $i $new
done

mkdir stripped

# Should do this on junc files
for i in $(ls GTEX*); do
	grep -v '^MT' $i > stripped/${i}
done

rm GTEX*

mv stripped/* .

# for i in $(ls *gz); do
# 	zgrep -v '^MT' $i > stripped/${i}.stripped
# 	bgzip stripped/${i}.stripped
# done

sbatch intronclustering.sh