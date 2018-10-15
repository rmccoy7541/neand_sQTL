# This program is to read and merge the neanderthal rsID data with the tissue sQTL data derived
# from using altrans on gTEX

library(data.table)
library(dplyr)

# args[1] will be the tissue file
# args[2] will be the neanderthal file
# args[3] maybe will be $SLURM_ARRAY_TASK_ID
args = commandArgs(trailingOnly=TRUE)

f1 <- fread(args[1])
f2 <- fread(args[2], fill=TRUE)
f2 <- f2[Num != 'PosID' | PosID != 'rsID']
# f1 <- f1[,-c(2,3)]
f2 <- f2[,-c("Num")]
colnames(f1)[colnames(f1)=="Variation"] <- "rsID"

merge_data <- merge(f2, f1, by = "rsID")

outfile <- paste(substr(args[1], 1, 3),"_", substr(args[2],1,5), "_output.txt", sep = "")
write.table(merge_data, file = outfile)

# merge_data <- merge(f2, f1[,c("Variation", "Direction", "Spearman_rho", "log10(p)")], rsID, by = "PosID")
