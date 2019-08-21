#!/bin/bash
#SBATCH --job-name=callJuncCluster
#SBATCH --time=1:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
# number of tasks (processes) per node
#SBATCH --ntasks-per-node=1
#######################################

# strips filenames of irrelevant info
for i in $(ls *junc); do
	new=$(echo $i | awk -F'.' '{print $1}')
	mv $i $new
done

mkdir stripped

# removing mitochondrial sequences
for i in $(ls GTEX*); do
	grep -v '^MT' $i > stripped/${i}
done


rm GTEX*

mv stripped/* .

rmdir stripped/
ls GTEX* | sort -V > juncfiles.txt