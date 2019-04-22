#!/bin/bash
#SBATCH --job-name=callJuncCluster
#SBATCH --time=1:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=1
#######################################

for i in $(ls *junc); do
	new=$(echo $i | awk -F'.' '{print $1}')
	mv $i $new
done

mkdir stripped

time {
# Should do this on junc files
for i in $(ls GTEX*); do
	grep -v '^MT' $i > stripped/${i}
done
}

rm GTEX*

mv stripped/* .

rmdir stripped/

# for i in $(ls *gz); do
# 	zgrep -v '^MT' $i > stripped/${i}.stripped
# 	bgzip stripped/${i}.stripped
# done

sbatch ${scripts}/../AWS/intronclustering.sh