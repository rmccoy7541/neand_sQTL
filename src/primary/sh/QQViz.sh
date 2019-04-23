#!/bin/bash

#SBATCH --partition=shared
#SBATCH --job-name=QQViz
#SBATCH --nodes=1
#SBATCH --time=1:0:0
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1

######################
# Begin work section #
######################

ml R
ml gcc

# this project's scripts dir
scripts=$(echo /home-1/aseyedi2@jhu.edu/work/aseyedi2/neand_sQTL/src/primary/)
# data dir
data=$(echo /home-1/aseyedi2@jhu.edu/work/aseyedi2/neand_sQTL/data/)


line=`sed "${SLURM_ARRAY_TASK_ID}q;d" GTExCovKey.csv`
full=$(echo $line | awk -F',' '{print $1}')
abb=$(echo $line | awk -F',' '{print $2}')
if grep "$abb" tissuesused.txt; then

   Rscript ${scripts}/R/QQPlot-Viz.R /home-1/aseyedi2@jhu.edu/work/aseyedi2/sQTL/$abb $abb/$abb.nominals.all.chunks.NE_only.txt.gz $abb/$abb.permutations_full.txt.gz ${data}/../analysis/SPRIME/sprime_calls.txt
   echo "$full is done with QQViz"
else
   echo "$full was not included in the analysis"
fi
