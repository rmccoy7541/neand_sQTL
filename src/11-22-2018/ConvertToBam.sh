# This converts the RNA-Seq GTEx data (.sra) to .bam files in preparation for using LeafCutter

for f in $PWD/*.sra; do ./sam-dump $f | samtools view -bS > $f.bam; done

