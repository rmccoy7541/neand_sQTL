library(data.table)
library(tidyverse)

file_path <- list.files(path = "~/Downloads/GTEx_Analysis_v8_eQTL_covariates/", 
                        pattern = ".txt",
                        full.names = TRUE)

read_tissue_covar <- function(file_path) {
  fread(file_path) %>%
    .[1:2,] %>%
    t(.) %>%
    .[-1,] %>%
    as.data.table(keep.rownames = TRUE) %>%
    setnames(., c("subject_id", "PC1", "PC2")) %>%
    .[, PC1 := as.numeric(PC1)] %>%
    .[, PC2 := as.numeric(PC2)] %>%
    return(.)
}

covar <- rbindlist(lapply(file_path, function(x) read_tissue_covar(x))) %>%
  .[!duplicated(subject_id)]

ggplot(data = covar, aes(x = PC1, y = PC2)) +
  geom_point()

covar[, cluster := kmeans(covar[, 2:3, with = FALSE], centers = 3)$cluster]

ggplot(data = covar, aes(x = PC1, y = PC2, color = factor(cluster))) +
  geom_point()

fwrite(covar[cluster != 3][, "subject_id"], file = "~/Downloads/exclude_samples.txt", 
       quote = FALSE, row.names = FALSE, col.names = FALSE)

