require("data.table")
require("R.utils")
args = commandArgs(trailingOnly=TRUE)
# args[1] is the leafcutter-generated phenotypes, args[2] is the tissue table
NE <- fread(paste0("zcat ", args[1]))
#NE <- fread("zcat Ne-sQTL_perind.counts.gz.qqnorm_chr16.gz.qtltools")
tistab <- fread(args[2])
#tistab <- fread("tissue_table.txt")
# below takes the SRR IDs found in NE column headers, matches them to those found in the
# tissue table, and then changes them the GTEX sample ID
ind <- match(names(NE), tistab$Run)
names(NE) <- tistab$submitted_subject_id[ind]

# I guess now what I want to do is find the tissue that corresponds to each sample, and write
# to file the phenotypes or whatever

## From GitHub
# If you split Sample_Name by body_site, you get a vector of Sample_Names corresponding to
# each body_site. Then you just need to intersect this with the names of NE for each body_site,
# and select the columns resulting from that intersection. The result is a named list of
# data tables. The names are the body_site values.
sites <- with(tistab, split(submitted_subject_id, body_site))

keep <- c('#Chr', 'start', 'end', 'ID')

tissues <- lapply(sites, function(x)
  NE[, .SD, .SDcols = c(keep, intersect(names(NE), x))])

sapply(names(tissues), function (x) write.table(tissues[[x]], file=paste0(paste0(tissues[[x]]$`#Chr`[1], "_"), x,".", "txt"), row.names=F, quote=FALSE, sep="\t"))