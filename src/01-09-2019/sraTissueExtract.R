library("data.table")
args = commandArgs(trailingOnly=TRUE)
# SraRunTable.txt is args[1]
sratabl <- fread(args[1])
# files <- fread("srafiles.txt", header=FALSE)
tissuesite <- subset(sratabl, select=c("Run", "Sample_Name", "body_site"))
write.table(tissuesite, file ="tissue_table", quote = F, sep= "\t", eol = "\r\n", row.names = F)