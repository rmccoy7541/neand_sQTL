library(data.table)
library(tidyverse)
library(ggplot2)
library(dplyr)

setwd(snakemake@input[["wd"]])
# countCounts
file_list <- list.files(pattern = "*.txt")

#read them all to file
for (file in file_list){
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- fread(file, header=TRUE, sep="\t")
    dataset$tissue <- gsub("^([^_]*_[^_]*)_.*$", "\\1", file)
  }
  # if the merged dataset does exist, append to it
  if (exists("dataset")){
    temp_dataset <-fread(file, header=TRUE, sep="\t")
    temp_dataset$tissue <- gsub("^([^_]*_[^_]*)_.*$", "\\1", file)
    # browser()
    if(nrow(temp_dataset) == 0 ){
      next();
    }
    dataset<-rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
}

dataset[is.na(dataset)] <- 0

# transcripts where the number of individuals with the HH not expressing the transcript is greater than zero and the
# number of individuals expressing the transcript with the same genotype is less than 10
loosened <- dataset[HH0 > 0 & `HH>0` < 10 & `HN>0` > 0][, .(n=.N), by=.(transcript_id)][order(n),]

df_uniq <- unique(loosened$transcript_id)

write.table(loosened,
            file = "results/loosenedRestrictions.txt",
            sep = "\t",
            row.names = F,
            quote = FALSE)
