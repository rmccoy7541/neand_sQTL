# Order sQTLs and transcripts
library(dplyr)
library(data.table)

filelist <- list.files(pattern = "*_filtered_iso.txt")

datalist <- lapply(filelist, FUN=fread, header = T)

for (i in 1:49){datalist[[i]]$tissue <- gsub("^([^_]*_[^_]*)_.*$", "\\1", filelist[i])}

datatab <- as.data.table(do.call(rbind, datalist))

datatab[order(variant_id, transcript_id, is_NL != 0)]
