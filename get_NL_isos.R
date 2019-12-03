library(tidyr)
library(data.table)

dt <- fread("Combined2_Ovary.txt")
dt1 <- dt[,17:197]
dt2 <- dt$Name
dt2 <- cbind(dt2, dt[,198:ncol(dt)])

dt1 <- as.data.table(dt1 %>% pivot_longer(-Name, names_to = "tissue_id", values_to = "counts"))

dt2 <- as.data.table(dt2 %>% pivot_longer(-Name, names_to = "individual", values_to = "is_NL"))

dt1$individual <- gsub("^([^.]*.[^.]*)..*$", "\\1", dt1$tissue_id)

final <- as.data.table(dplyr::full_join(dt1, dt2, by = c("Name", "individual")))

final <- final[is_NL == 1] 