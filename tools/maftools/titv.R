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
   "--output_pdf", "-o",
   required="True",
   help="Output directory to store gene plots"
   )

args <- parser$parse_args();

###

laml = read.maf(maf = args$input_maf, removeSilent = T, useAll = F)
pdf(args$output_pdf)
titv(maf = laml, useSyn=T)
dev.off()
