library(data.table)
library(tidyverse)

args = commandArgs(trailingOnly=TRUE)

dt <- fread(args[1], header = T)
colnames(dt) <- c("variant_id", "transcript_id", "HH", "HNL", "NLNL")

dt <- dt %>%
  mutate_at(vars(HH:NLNL), list(ratio = ~str_split(., ";") )) %>%
  mutate_at(vars(HH_ratio:NLNL_ratio), list(~map_dbl(., ~as.numeric(.x) %>%
  {.[[1]]/.[[2]]})))

dt$diff <- (dt$NLNL_ratio - dt$HH_ratio)

tissue_name <- gsub("^([^_]*_[^_]*)_.*$", "\\1", args[1])

write.table(dt,
            file = paste0(tissue_name, "_filteredTable_horizontal_difference.txt"),
            sep = "\t",
            row.names = F,
            quote = FALSE)