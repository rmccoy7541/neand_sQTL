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

plot2D(termCluster, 1)