library(data.table)
library(tidyverse)

dt <- fread("Muscle_Skeletal_filteredTable_horizontal.txt", header = T)
colnames(dt) <- c("variant_id", "transcript_id", "HH", "HNL", "NLNL")

dt <- dt %>%
  mutate_at(vars(HH:NLNL), list(new = ~str_split(., ";") )) %>%
  mutate_at(vars(HH_new:NLNL_new), list(~map_dbl(., ~as.numeric(.x) %>%
  {.[[1]]/.[[2]]})))
