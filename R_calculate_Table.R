library(data.table)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

tissue_name <- args[1]

dt <- fread(args[1], header = T)
colnames(dt) <- c("variant_id", "transcript_id", "HH", "HNL", "NLNL")

dt <- dt %>%
  na.omit() %>%
  mutate_at(vars(HH:NLNL), list(ratio = ~str_split(., ";") )) %>%
  mutate_at(vars(HH_ratio:NLNL_ratio), list(~map_dbl(., ~as.numeric(.x) %>%
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