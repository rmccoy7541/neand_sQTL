library("data.table")
args = commandArgs(trailingOnly=TRUE)
# SraRunTable.txt is args[1]
sratabl <- fread(args[1])
# current directory is args [2]
tissuesite <- subset(sratabl, select=c("Run", "Sample_Name", "body_site", "submitted_subject_id"))
write.table(tissuesite, file=paste(args[2], "tissue_table.txt"), quote = F, sep= "\t", eol = "\r\n", row.names = F)