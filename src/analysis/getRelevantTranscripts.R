library(data.table)
library(tidyverse)
library(ggplot2)

file_list <- list.files(pattern = "*.txt")

# *horizontal
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

theme_set(theme_bw())  # pre-set the bw theme.

otherdt <- vector(mode = "numeric")

## make it so that the cutoff value includes the sum of the unique transcripts for cut off values below it
for(i in 0:10) {
  if (i == 0){
  otherdt[i + 1] <- sum(dataset[HH0 > 0 & `HH>0` == i & `HN>0` > 0 & HN0 == 0][, .(n=.N), by=.(transcript_id)][order(n),]$n)
  }
  else {
  otherdt[i + 1] <- sum(sum(dataset[HH0 > 0 & `HH>0` == i & `HN>0` > 0 & HN0 == 0][, .(n=.N), by=.(transcript_id)][order(n),]$n), otherdt[i])
  }
}

# transcripts where the number of individuals with the HH not expressing the transcript is greater than zero and the
# number of individuals expressing the transcript with the same genotype is less than 10
loosened <- dataset[HH0 > 0 & `HH>0`< 10][, .(n=.N), by=.(transcript_id)][order(n),]

write.table(loosened,
            file = "loosenedRestrictions.txt",
            sep = "\t",
            row.names = F,
            quote = FALSE)
