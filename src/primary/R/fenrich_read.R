filenames <- list.files(".", pattern="*enrichment.QTL.txt", full.names=TRUE)
D=do.call(rbind, lapply(filenames, function(x) read.table(x, header=FALSE, stringsAsFactors=FALSE, col.names=c("observed.QTLs", "total.QTLs", "expected.QTLs", "stdDev.QTLs"))))

toKeep <- seq(1, nrow(D), 2)
D <- D[toKeep,]

for(i in 1:47){
  D[i,5] <-   gsub("/", "", strsplit(filenames[i], "\\.")[[1]][2])
}
rownames(D) <- D$V5
D <- D[,-5]

ft.res <- apply(D, 1, function(x){
    t1 <- fisher.test(matrix(x, nrow = 2))
    data.frame(p_value = t1$p.value, odds_ratio = t1$estimate)
})

D <- cbind(D, do.call(rbind, ft.res))

write.table(D, file = "fenrich_table.txt", quote = F, sep= "\t", eol = "\r\n", row.names = F)