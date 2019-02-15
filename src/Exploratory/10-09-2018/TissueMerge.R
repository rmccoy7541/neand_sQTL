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
# not sure what I did here... I need someone's explanation
f2 <- f2[Num != 'PosID' | PosID != 'rsID']
<<<<<<< HEAD
<<<<<<< HEAD
# Remove ensembl IDs - removed it on 10/15/18 to retain them
# f1 <- f1[,-c(2,3)]
=======
=======
>>>>>>> 97d3a0806c2cc02fbc7bdf11e0e09c70eeb4016c
# making the names unique because the two ensembl ID columns are both named 'Link'
names(f1) <- make.unique(names(f1))
#Next line removes ensembl IDs - removed that on 10/15/2018
#f1 <- f1[,-c(2,3)]
<<<<<<< HEAD
>>>>>>> 97d3a0806c2cc02fbc7bdf11e0e09c70eeb4016c
=======
>>>>>>> 97d3a0806c2cc02fbc7bdf11e0e09c70eeb4016c
f2 <- f2[,-c("Num")]
colnames(f1)[colnames(f1)=="Variation"] <- "rsID"

merge_data <- merge(f2, f1, by = "rsID")

outfile <- paste(substr(args[1], 9, 12),"_", substr(args[2],9,11), "_output.txt", sep = "")
write.table(merge_data, file = outfile)

# merge_data <- merge(f2, f1[,c("Variation", "Direction", "Spearman_rho", "log10(p)")], rsID, by = "PosID")
