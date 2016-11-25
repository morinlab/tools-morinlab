maf <- read.delim("./MAF", header=T)
maf1 <- maf[,c("Tumor_Sample_Barcode", "Gene", "SIFT", "PolyPhen")]
write.table(maf1, quote = FALSE, row.names = FALSE) 
