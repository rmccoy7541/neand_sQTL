# 12-10-2018

## filter_bam.sh
This script filters from the `.bam` files all unplaced contigs and basically anything that isn't found on chromosomes 1-22, sex or mitochondrial. As of today, it is unclear if the 'M' chromosomes can be handled by LeafCutter. We shall see. 

## bam2junccall.sh
This script sets up the SLURM job array in which LeafCutter's `bam2junc.sh` script is called and used on each filtered `.bam` file. Remember to set full path of LeafCutter's directory in `bam2junc.sh`.
