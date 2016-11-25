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
   "--top", "-t",
   required="True",
   help="The Top [val] genes to include in plot"
   )

parser$add_argument(
    "--gene_mask_list", "-sl",nargs='?',
     help="Input gene blacklist separated by newline (optional)", default=0
			  );
						

parser$add_argument(
   "--output_pdf", "-o",
   required="True",
   help="Output directory to store gene plots"
   )

args <- parser$parse_args();

###

if(args$gene_mask_list != 0){
all_genes <- read.table(args$gene_mask_list, stringsAsFactors=FALSE)[,1]
laml = read.maf(maf = args$input_maf, removeSilent = T, useAll = T)
pdf(args$output_pdf)
oncoplot(maf = laml, top=args$top,genesToIgnore=all_genes)
dev.off()
}else{
laml = read.maf(maf = args$input_maf, removeSilent = T, useAll = T)
pdf(args$output_pdf)
oncoplot(maf = laml, top=args$top)
dev.off()
}