library(tidyr)
library(data.table)

dt <- fread("combined2_Ovary.txt")

list_dt <- split.default(dt, nchar(names(dt)) > 10)

xcrips <- cbind(list_dt[[1]][,"variant_id"], list_dt[[1]][,"Name"], list_dt[[2]][,5:ncol(list_dt[[2]])])

nl_iso <- cbind(list_dt[[1]][,"variant_id"], list_dt[[1]][,"Name"], list_dt[[1]][,13:ncol(list_dt[[1]])])

#add variant id
xcrips <- as.data.table(xcrips %>% pivot_longer(-c(Name, variant_id), names_to = "tissue_id", values_to = "counts"))

#add variant id
nl_iso <- as.data.table(nl_iso %>% pivot_longer(-c(Name, variant_id), names_to = "individual", values_to = "is_NL"))

xcrips$individual <- gsub("^([^.]*.[^.]*)..*$", "\\1", xcrips$tissue_id)

final <- as.data.table(dplyr::full_join(xcrips, nl_iso, by = c("Name", "individual")))
