args <- commandArgs(trailingOnly = TRUE)

files <- read.table(args[1], stringsAsFactors=F)

my.read.table <- function(file) {
    read.table(
        file,
        sep="\t",
        stringsAsFactors=T,
        header=as.logical(args[4])
        )
    }

data <- NA;
if (length(files) == 1) {
    data <- my.read.table(files[1])
    }
else {
    data_list <- apply(files, 1, my.read.table)
    data <- Reduce(function(x,y) { rbind(x,y) }, data_list)
    }

contig_order <- read.table(args[2], stringsAsFactors=F)[,1]
data[,1] <- factor(data[,1], levels=contig_order)
data <- data[order(data[,1]),]
data[,1] <- as.character(data[,1])

write.table(
    x=data,
    file=args[3],
    sep ="\t",
    quote=F,
    row.names=F,
    col.names=as.logical(args[4])
    )
