require("data.table")
requirse("R.utils")
args = commandArgs(trailingOnly=TRUE)
# args[1] is the leafcutter-generated phenotypes, args[2] is the tissue table
NE <- fread("NE_sQTL_perind.counts.gz.qqnorm_chr13")

tistab <- fread("tissue_table.txt")

# below takes the SRR IDs found in NE column headers, matches them to those found in the
# tissue table, and then changes them the GTEX sample ID
ind <- match(names(NE), tistab$Run)
names(NE) <- tistab$Sample_Name[ind]

# I guess now what I want to do is find the tissue that corresponds to each sample, and write
# to file the phenotypes or whatever

## From GitHub (https://stackoverflow.com/questions/54189236/how-do-i-split-a-table-into-several-new-tables-based-on-whether-the-column-heade/54189469#54189469)

# If you split Sample_Name by body_site, you get a vector of Sample_Names corresponding to 
# each body_site. Then you just need to intersect this with the names of NE for each body_site, 
# and select the columns resulting from that intersection. The result is a named list of 
# data tables. The names are the body_site values.
sites <- with(tistab, split(Sample_Name, body_site))

keep <- c('#Chr', 'start', 'end', 'ID')

lapply(sites, function(x) 
  NE[, .SD, .SDcols = c(keep, intersect(names(NE), x))])