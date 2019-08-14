require("data.table")
require("R.utils")
args = commandArgs(trailingOnly=TRUE)
#args[1] is the leafcutter-generated phenotypes, args[2] is the tissue table
NE <- fread(paste0("zcat ", args[1]))

setnames(NE, c('ID', '.', '+'), c('PID', 'GID', 'Strand'))

tistab <- fread(args[2])

# below takes the SRR IDs found in NE column headers, matches them to those found in the
# tissue table, and then changes them the GTEX sample ID
ind <- match(names(NE), tistab$Run)
names(NE) <- tistab$Sample_Name[ind]

# I guess now what I want to do is find the tissue that corresponds to each sample, and write
# to file the phenotypes or whatever

## From GitHub
# If you split Sample_Name by body_site, you get a vector of Sample_Names corresponding to
# each body_site. Then you just need to intersect this with the names of NE for each body_site,
# and select the columns resulting from that intersection. The result is a named list of
# data tables. The names are the body_site values.
sites <- with(tistab, split(Sample_Name, body_site))

keep <- c('#Chr', 'start', 'end', 'PID', 'GID', 'Strand')

tissues <- lapply(sites, function(x)
  NE[, .SD, .SDcols = c(keep, intersect(names(NE), x))])

## keep only the tables which have more than 6 columns
keep_index <- unlist(lapply(names(tissues), function(x) ncol(tissues[[eval(quote(x))]]) > 6))
tissues_kept <- tissues[keep_index]

# change the headers back to subject id
for (i in 1:length(tissues_kept)){
  ind <- match(names(tissues_kept[[i]]), tistab$Sample_Name)
  names(tissues_kept[[i]]) <- tistab$submitted_subject_id[ind]
}


sapply(names(tissues_kept), function (x) write.table(tissues_kept[[x]], file=paste0(paste0(tissues_kept[[x]]$`#Chr`[1], "_"), x,".", "txt"), row.names=F, quote=FALSE, sep="\t"))
