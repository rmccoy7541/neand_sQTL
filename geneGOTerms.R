library(data.table)
library(mygene)
library(ggplot2)
library(RDAVIDWebService)
library(biomaRt)

test <- fread("results/loosenedRestrictionsGenes.txt")

ensembl_id <- gsub("\\..*","",test$gene_id)

david <- DAVIDWebService$new(email="aseyedi2@jhu.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

result <- addList(david, ensembl_id, idType = "ENSEMBL_GENE_ID", listName = "NL_genes", listType = "Gene")

setAnnotationCategories(david, c("GOTERM_BP_ALL"))

termCluster<-getClusterReport(david, type="Term")

# List of genes in cluster 1
# 
ensembl_genes <- unlist(strsplit(members(termCluster)[[1]][[6]][6], split=", "))

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

#list of genes with highest DAVID "enrichment" value and 
G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","uniprot_gn_symbol", "description"),values=ensembl_genes,mart= mart)

write.table(G_list,
            file = paste0("results/NL_Iso_TopGenes_List.txt"),
            sep = "\t",
            row.names = F,
            quote = FALSE)