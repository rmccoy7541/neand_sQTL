library("data.table")
# args = commandArgs(trailingOnly=TRUE)
# SraRunTable.txt is args[1], tissue key is args[2]
sratabl <- fread(snakemake@input[[1]])
tiskey <- fread(snakemake@input[[2]])

tissuesite <- subset(sratabl, select=c("Run", "Sample_Name", "body_site", "submitted_subject_id"))

# get the tissue sites for each corresonding sra file
tissuesite$tmp <- tiskey$Key[ match(tissuesite$body_site, tiskey$Tissue) ]
tissuesite$tmp <- ifelse(is.na(tissuesite$tmp), tissuesite$body_site, tissuesite$tmp)

tissuesite$body_site <- tissuesite$tmp
tissuesite$tmp <- NULL

write.table(tissuesite, file=gsub(" ", "", paste("tissue_table.txt")), quote = F, sep= "\t", eol = "\r\n", row.names = F)