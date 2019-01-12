library("data.table")
args = commandArgs(trailingOnly=TRUE)
# SraRunTable.txt is args[1]
# Current directory is args[2]
sratabl <- fread(args[1])
# files <- fread("srafiles.txt", header=FALSE)
tissuesite <- subset(sratabl, select=c("Run", "Sample_Name", "body_site"))
write.table(tissuesite, file=paste(args[2],"tissue_table.txt"), quote = F, sep= "\t", eol = "\r\n", row.names = F)