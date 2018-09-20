# Going to write an R script that joins the rsID data with
# the Neanderthal data on the basis of the position ID, which
# is the first column for both

library(data.table)

args = commandArgs(trailingOnly=TRUE)

NEdata <- fread(args[1], header=F)

names(NEdata) <- c("PosID", "aAllele", "derAllele")

rsID <- fread(args[2], header=F)

names(rsID) <- c("PosID", "rsID")

merge_data <- merge(NEdata[, c("PosID")], rsID, by = "PosID")

outfile <- paste(args[4],"_", args[3], "_output.txt", sep = "")
write.table(merge_data, file = outfile)