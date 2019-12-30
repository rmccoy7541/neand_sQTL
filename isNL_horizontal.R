library(tidyverse)
library(data.table)

#args[1] is filtered table,
args = commandArgs(trailingOnly=TRUE)

dt <- fread(args[1]) %>%
  mutate(counts = counts %>% as.character(),
  nrows = nrows %>% as.character()) %>% 
  unite("result",counts,nrows,sep = ";") %>% 
  pivot_wider(names_from = is_NL,values_from = result) %>%
  as.data.table()

tissue_name <- gsub("^([^_]*_[^_]*)_.*$", "\\1", args[1])

write.table(dt,
            file = paste0(tissue_name, "_filteredTable_horizontal.txt"),
            sep = "\t",
            row.names = F,
            quote = FALSE)