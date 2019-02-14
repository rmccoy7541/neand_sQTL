library(data.table)
setwd('D:/Dropbox/Dropbox/Grad school/McCoy Lab/matches')
EAS <- fread("EASmatch.txt", header = F)
EAS$V1 <- NULL