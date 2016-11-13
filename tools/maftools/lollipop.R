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
    "--gene_list", "-gl",
    required="True",
    help="Input gene list with separated by newline"
    );

parser$add_argument(
   "--output_directory", "-o",
   required="True",
   help="Output directory to store gene plots"
   )

args <- parser$parse_args();

###

all_genes <- read.table(args$gene_list, stringsAsFactors=FALSE)[,1]
laml = read.maf(maf = args$input_maf, removeSilent = T, useAll = T)
count = 1
for (gene in all_genes) {
    output_name = paste(args$output_directory, paste("samp", count, ".pdf", sep=""), sep="/")
    pdf(output_name)
    lollipopPlot(maf = laml, gene = gene, AACol = 'amino_acid_change', labelPos = 'all')
    dev.off()
    count <- count + 1
    }
