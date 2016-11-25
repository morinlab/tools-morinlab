maf <- read.delim("./MAF", header=TRUE, stringsAsFactors=FALSE, na.strings=c(""," ","NA"))
colnames(maf) = sapply(colnames(maf), tolower)
maf$sift[maf$variant_classification=="Frame_Shift_Del"] <- 0
maf$sift[maf$variant_classification=="Frame_Shift_Ins"] <- 0
maf$sift[maf$variant_classification=="Nonsense_Mutation"] <- 0
maf$sift[maf$variant_classification=="Translation_Start_Site"] <- 0
maf$sift[maf$variant_classification=="Splice_Site"] <- 0
maf$sift[maf$variant_classification=="Silent"] <- 1
maf$polyphen[maf$variant_classification=="Frame_Shift_Del"] <- 1
maf$polyphen[maf$variant_classification=="Frame_Shift_Ins"] <- 1
maf$polyphen[maf$variant_classification=="Nonsense_Mutation"] <- 1
maf$polyphen[maf$variant_classification=="Translation_Start_Site"] <- 1
maf$polyphen[maf$variant_classification=="Splice_Site"] <- 1
maf$polyphen[maf$variant_classification=="Silent"] <- 0
maf1 <- maf[,c("tumor_sample_barcode", "hugo_symbol", "sift", "polyphen")]
write.table(maf1, sep = '\t', quote = FALSE, row.names = FALSE)