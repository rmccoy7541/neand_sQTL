library(data.table)
library(mygene)
library(ggplot2)
library(RDAVIDWebService)
library(biomaRt)

test <- fread("results/loosenedRestrictionsGenes.txt")

#hg38
# ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
# 
# ensemblIDs1 <- getBM(attributes='ensembl_gene_id', 
#                 filters = 'hgnc_symbol', 
#                 values = test$gene_name, 
#                 mart = ensembl)

david <- DAVIDWebService$new(email="aseyedi2@jhu.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")

ei <- fread("entrez_id_conversion.txt")

result <- addList(david, ei$To, idType = "ENTREZ_GENE_ID", listName = "NL_genes", listType = "Gene")

setAnnotationCategories(david, c("GOTERM_BP_ALL"))

termCluster<-getClusterReport(david, type="Term")

plot2D(termCluster, 1)
# this is basically unapproachably challenging for me because 
# a) the number of items in res[[1]] (res>listData>query) is 392 even though the number of gene_names in test is 384 AND the number of unique items in both is 144. 
# I don't know what is going on there.
# b) I need, within res, each list element in go.BP to be matched to the corresponding gene symbol in query, and then those queries to be matched to the corresponding 
# transcripts in test. And I don't know how to do that. I don't even know how to really articulate how to do this to someone who doesn't have the files I'm talking about.
# So this might be another one for Rajiv. Although I could probably figure out how to do this if I just tried to psychically bend the spoon hard enough, I think I would
# rather tackle easier things first and then address this later, possibly with Rajiv.



