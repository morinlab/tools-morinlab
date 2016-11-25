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
     type="integer",
    help="Minimum number of mutations for each included gene"
    )

parser$add_argument(
    "--gene_mask_list", "-sl",nargs='?',
	    help="Input gene blacklist separated by newline (optional)", default=0
	        );

parser$add_argument(
    "--output_pdf", "-o",
    required="True",
    help="Output PDF"
    )

args <- parser$parse_args();

###



if(args$gene_mask_list != 0){

	all_genes <- read.table(args$gene_mask_list, stringsAsFactors=FALSE)[,1]


	laml = read.maf(maf = args$input_maf, removeSilent = T, useAll = T)
	pdf(args$output_pdf)
	geneCloud(input = laml, minMut = args$min_mutations,genesToIgnore=all_genes)
	dev.off()
} else{

  laml = read.maf(maf = args$input_maf, removeSilent = T, useAll = T)
  pdf(args$output_pdf)
  geneCloud(input = laml, minMut = args$min_mutations)
  dev.off()
}