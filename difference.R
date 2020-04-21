library(data.table)
library(tidyverse)

# args = commandArgs(trailingOnly=TRUE)

# Horizontalize
dt <- fread("Whole_Blood_NL_isos.txt") %>%
  mutate(counts = counts %>% as.character(),
         nrows = nrows %>% as.character()) %>%
  unite("result",counts,nrows,sep = ";") %>%
  pivot_wider(names_from = is_NL,values_from = result) %>%
  as.data.table()

# get tissue name to put in a column
tissue_name <- gsub("^([^_]*_[^_]*)_.*$", "\\1", "Whole_Blood_NL_isos.txt")

colnames(dt) <- c("variant_id", "transcript_id", "individual", "HH", "HN", "NN")

dt <- dt %>%
  na.omit() %>%
  mutate_at(vars(HH:NN), list(ratio = ~str_split(., ";") )) %>%
  mutate_at(vars(HH_ratio:NN_ratio), list(~map_dbl(., ~as.numeric(.x) %>%
                                                       {.[[1]]/.[[2]]}))) %>%
  as.data.table()

dt$diff <- (dt$NLNL_ratio - dt$HH_ratio)

setorder(dt, -diff)

tissue_name <- gsub("^([^_]*_[^_]*)_.*$", "\\1", tissue_name)

write.table(dt,
            file = paste0(tissue_name, "_horizontal_difference.txt"),
            sep = "\t",
            row.names = F,
            quote = FALSE)