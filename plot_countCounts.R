library(data.table)
library(tidyverse)
library(ggplot2)

file_list <- list.files(pattern = "_countCounts.txt", full.names = TRUE)
f_dowle2 = function(DT) {
  for (i in names(DT))
    DT[is.na(get(i)), (i):=0]
}
read_tissue_counts <- function(filename) {
  dt_tissue <- fread(filename, header = TRUE) 
  f_dowle2(dt_tissue)
  tissue_name <- gsub("^([^_]*_[^_]*)_.*$", "\\1", filename)
  tissue_name <- gsub("/Users/rajivmccoy/Downloads/", "", tissue_name)
  dt_tissue[, tissue := tissue_name]
  return(dt_tissue[!duplicated(transcript_id)])
}
dt <- do.call(rbind, lapply(file_list, function(x) read_tissue_counts(x)))

ggplot(dt, aes(x = `HH>0`, y = `HN>0` + `NN>0`)) +
  geom_bin2d() + 
  scale_fill_gradient(name = "count", trans = "log10") +
  facet_wrap(~ tissue)

ggplot(dt, aes(x = `HH>0`, y = `HN>0` + `NN>0`)) +
  geom_bin2d(bins = 100) + 
  scale_fill_gradient(name = "count", trans = "log10")



# file_list <- list.files(pattern = "*.txt")
# 
# # # *horizontal
# # for (file in file_list){
# #   # if the merged dataset doesn't exist, create it
# #   if (!exists("dataset")){
# #     dataset <- fread(file, header=TRUE, sep="\t")
# #     dataset$tissue <- gsub("^([^_]*_[^_]*)_.*$", "\\1", file)
# #   }   
# #   # if the merged dataset does exist, append to it
# #   if (exists("dataset")){
# #     temp_dataset <-fread(file, header=TRUE, sep="\t")
# #     temp_dataset$tissue <- gsub("^([^_]*_[^_]*)_.*$", "\\1", file)
# #     # browser()
# #     if(nrow(temp_dataset) == 0 ){
# #       next();
# #     }
# #     dataset<-rbind(dataset, temp_dataset)
# #     rm(temp_dataset)
# #   }
# # }
# 
# dataset <- fread("Whole_Blood_countCounts.txt", header=TRUE, sep="\t")
# 
# dataset[is.na(dataset)] <- 0
# 
# theme_set(theme_bw())  # pre-set the bw theme.
# 
# otherdt <- vector(mode = "numeric")
# 
# 
# ## make it so that the cutoff value includes the sum of the unique transcripts for cut off values below it
# for(i in 0:10) {
#   otherdt[i + 1] <- sum(dataset[HH0 > 0 & `HH>0` == i & `HN>0` > 0][, .(n=.N), by=.(transcript_id)][order(n),]$n)
# }
# 
# sum(dataset[HH0 > 0 & `HH>0` == 0][, .(n=.N), by=.(transcript_id)][order(n),]$n)
# 
# g <- ggplot(dataset, aes(transcript_id, HH0))
# 
# # Scatterplot
# g + geom_point() + 
#   geom_smooth(method="lm", se=F) +
#   labs(subtitle="# of uniq transcripts for each transcript x sQTL pair with ascending HH counts", 
#        y="xcript_id", 
#        x="HH0", 
#        title="Scatterplot with overlapping points")