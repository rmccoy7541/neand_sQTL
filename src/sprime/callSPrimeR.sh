#!/bin/bash
#SBATCH --job-name=SPrimeR
#SBATCH --time=4:00:00
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8

Rscript SPrime.R results.score $splice_ai