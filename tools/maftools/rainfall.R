require(maftools);
library(argparse);

###

parser <- ArgumentParser(description="Create a Gene Lollipop using Maftools");

parser$add_argument(
    "--input_maf", "-maf",
    required="True",
    help="Input Variants in MAF format"
    );

parser$add_argument(
    "--sample_list", "-sl",
    required="True",
    help="Input sample list separated by newline"
    );

parser$add_argument(
   "--output_directory", "-o",
   required="True",
   help="Output directory to store rainfall plots"
   )

args <- parser$parse_args();

###

all_samples <- read.table(args$sample_list, stringsAsFactors=FALSE)[,1]
laml = read.maf(maf = args$input_maf, removeSilent = T, useAll = T)
count = 1
for (sample in all_samples) {
    output_name = paste(args$output_directory, paste("samp", count, ".pdf", sep=""), sep="/")
    pdf(output_name)
    rainfallPlot(maf = laml, tsb = sample)
    dev.off()
    count <- count + 1
    }
