library(data.table)
## I need to write this script to reorder the column headers which are now apparently out of wack
## I just need to shift them over one
# local file directory
#setwd("D:/")
#setwd("Dropbox/Dropbox/GitHub/neanderthal-sqtl/data/10-15-2018/NE_altrans/EUR_altrans/")

# Each file is named into a character vector list "filelist", and the new names for the column headers
# are then assigned to a variabled. lapply is used to make an actual list using 
filelist <- list.files(pattern = ".*.txt")
my_names <- c("RowID", "rsID", "PosID", "Link", "Link.1","Direction", "Spearman_rho", "-log10(p)")
datalist <- lapply(filelist, fread, fill = TRUE, col.names = my_names)
# names(datalist) <- filelist
for(dt in datalist) {
  write.table(dt, file = paste("/out/", dt))
}
