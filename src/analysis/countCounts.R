library(data.table)
library(tidyverse)

#args[1] is the final iso table

dt <- fread(snakemake@input[["finalIso"]], header = T)

tissue_name <- gsub("^([^_]*_[^_]*)_.*$", "\\1", args[2])

dt$nrows <- NULL

dt <- dt %>% 
  group_by(variant_id,transcript_id,is_NL) %>% 
  mutate(grp = case_when(counts == 0 ~ "0", counts > 0 ~ ">0")) %>% 
  count(grp) %>%
  ungroup() %>%
  mutate(
    is_NL = case_when(
      is_NL == 0 ~ "HH",
      is_NL == 1 ~ "HN",
      is_NL == 2 ~ "NN",
    )
  ) %>% 
  unite(comb,is_NL,grp,sep="") %>%
  pivot_wider(names_from = comb, values_from = n)


write.table(dt,
            file = paste0("results/countCounts/", tissue_name, "_countCounts.txt"),
            sep = "\t",
            row.names = F,
            quote = FALSE)