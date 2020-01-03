library(data.table)
library(dplyr)

tissuename <- snakemake@input["tisnames"]

iso <- fread(snakemake@input["isos"])

iso <- setDT(iso)[with(iso, (is_NL == 0 & counts/nrows < 50) | (is_NL %in% c(1,2) & counts/nrows > 50)),][, triplet := .N, by = .(variant_id, transcript_id)][triplet == 3, ][, triplet := NULL]

write.table(iso,
            file = paste0(tissuename, "_filtered_iso.txt"),
            sep = "\t",
            row.names = F,
            quote = FALSE)