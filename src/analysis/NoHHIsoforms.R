library(data.table)
library(tidyverse)

file_list <- list.files(pattern = "*al.txt")

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

ascendingGT <- dataset[grepl("^0", `0`) & !grepl('^0', `1`) & !grepl('^0', `2`)]

write.table(ascendingGT,
            file = paste0("NoHHIsoforms.txt"),
            sep = "\t",
            row.names = F,
            quote = FALSE)
