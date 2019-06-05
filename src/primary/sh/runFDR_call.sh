ml R
ml gcc

for i in $(ls *permutations.txt); do
	

	
/scratch/groups/rmccoy22/progs/QTLtools/script/runFDR_cis.R THYROID_permutations.txt 0.10 THYROID.results