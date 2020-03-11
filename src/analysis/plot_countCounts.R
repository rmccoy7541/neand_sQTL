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

# dataset <- fread("Whole_Blood_countCounts.txt", header=TRUE, sep="\t")
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
# otherdt <- sum(dataset[HH0 > 0 & `HH>0` == 0 & NN0 == 0 & `HN>0` >0 & `NN>0` > 0 & NN0 == 0][, .(n=.N), by=.(transcript_id)][order(n),]$n)
otherdt

sum(dataset[HH0 > 0 & `HH>0` == 0][, .(n=.N), by=.(transcript_id)][order(n),]$n)

g <- ggplot(dataset, aes(transcript_id, HH0))

# Scatterplot
g + geom_point() + 
  geom_smooth(method="lm", se=F) +
  labs(subtitle="# of uniq transcripts for each transcript x sQTL pair with ascending HH counts", 
       y="xcript_id", 
       x="HH0", 
       title="Scatterplot with overlapping points")