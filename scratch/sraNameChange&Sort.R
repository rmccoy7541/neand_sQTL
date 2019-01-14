library("data.table")
args = commandArgs(trailingOnly=TRUE)
# args[1] is the leafcutter-generated PCs, args[2] is the tissue table
NEPCs <- fread("NE_sQTL_perind.counts.gz.PCs")
tistab <- fread("tissue_table.txt")
# what do i want to do? rename each column in NEPCs based on its full GTEx ID. I think that's it. 