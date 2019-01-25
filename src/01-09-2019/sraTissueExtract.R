library("data.table")
args = commandArgs(trailingOnly=TRUE)
# SraRunTable.txt is args[1]
sratabl <- fread(args[1])
#sratabl <- fread("SraRunTable.txt")
#tiskey <- fread("GTExTissueKey.csv")
tiskey <- fread(args[2])

tissuesite <- subset(sratabl, select=c("Run", "Sample_Name", "body_site", "submitted_subject_id"))

tissuesite$tmp <- tiskey$Key[ match(tissuesite$body_site, tiskey$Tissue) ]
tissuesite$tmp <- ifelse(is.na(tissuesite$tmp), tissuesite$body_site, tissuesite$tmp)

tissuesite$body_site <- tissuesite$tmp
tissuesite$tmp <- NULL
#new <- sratabl
#new[] <- tiskey$Key[match(unlist(sratabl), tiskey$Tissue)]

write.table(tissuesite, file=gsub(" ", "", paste("tissue_table.txt")), quote = F, sep= "\t", eol = "\r\n", row.names = F)