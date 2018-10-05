# This program is to read and merge the neanderthal rsID data with the tissue sQTL data derived
# from using altrans on gTEX

library(data.table)
library(dplyr)

f1 <- fread("AdiposeSubcutaneous.Altrans.bestPerLink")
f2 <- fread("file:///D:/neanderthal_foolinaround/matches/EASmatch.txt", fill=TRUE)
f2 <- f2[Num != 'PosID' | PosID != 'rsID']
f1 <- f1[,-c(2,3)]
f2 <- f2[,-c("Num")]
colnames(f1)[colnames(f1)=="Variation"] <- "rsID"

merge_data <- merge(f2, f1, by = "rsID")

# merge_data <- merge(f2, f1[,c("Variation", "Direction", "Spearman_rho", "log10(p)")], rsID, by = "PosID")
