# Bash

Contained in this folder are Bash scripts, which were used to perform various, undemanding tasks. I will now describe them one-by-one.

## callNEjoinR.sh

This script invokes NEjoin.R, which matched the rsID data with the neanderthal posID data. It uses the various tools provided by SLURM to match each chromosome to the rsID
data independently, and passes command line arguments to the R script.

## charScript.sh

Very inefficently extracted the relevant data from the 1kg files.

## ChrCat.sh

Contatenated the chromosome number with the posID of each of the chromosome files (columns 1, 2 and 3, respectively).

## NEjoin.sh

The original shell version of the R script that I ended up using instead.

## NEsort.sh

Sorted the neanderthal data for the `join` call that I had originally used to match the neanderthal data with the rsIDs.

## sort.sh

Sorted rsID data according to first column, in preparation for calling `NEjoin.sh`.