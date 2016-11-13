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
    "--min_mutations", "-m",
    required="True",
    help="Minimum number of mutations for each included gene"
    )

parser$add_argument(
    "--output_pdf", "-o",
    required="True",
    help="Output PDF"
    )

args <- parser$parse_args();

###

laml = read.maf(maf = args$input_maf, removeSilent = T, useAll = T)
pdf(args$output_pdf)
geneCloud(input = laml, minMut = args$min_mutations)
dev.off()
