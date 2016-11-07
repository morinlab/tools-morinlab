require(maftools);
library(argparse);
library(NMF);

###

parser <- ArgumentParser(description="Create a Gene Lollipop using Maftools");

parser$add_argument(
    "--input_maf", "-maf",
    required="True",
    help="Input Variants in MAF format"
    );

parser$add_argument(
   "--genome", "-g",
   required="True",
   help="Reference Genome"
   )

parser$add_argument(
   "--output_pdf", "-o",
   required="True",
   help="Output directory to store gene plots"
   )

args <- parser$parse_args();

###

laml = read.maf(maf = args$input_maf, removeSilent = T, useAll = F)
tnm = trinucleotideMatrix(maf = laml, ref_genome = args$genome, prefix = '', add = T, useSyn = T)
sig = extractSignatures(mat = tnm)
pdf(args$output_pdf)
plotSignatures(nmfRes = sig)
dev.off()
