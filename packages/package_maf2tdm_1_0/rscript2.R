maf <- read.delim("./output", header=T, na.strings=c(""," ","NA"))
colnames(maf) = sapply(colnames(maf), tolower)
final <- maf[!(is.na(maf$sift)) | !(is.na(maf$polyphen)),]
write.table(final, sep = '\t', quote = FALSE, row.names = FALSE, col.names = TRUE)