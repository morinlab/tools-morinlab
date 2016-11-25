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
    "--gene_blacklist", "-gl",
    help="Input gene list with separated by newline"
    );

parser$add_argument(
   "--min_mut", "-mm",
      default=5,
         help="Minimum number of mutations seen in the gene for it to be included in the calculation");

parser$add_argument(
   "--aacol", "-ac",
   help="Optionally provide the name of the column that contains the amino acid annotation in your MAF file");

parser$add_argument(
   "--output_detail", "-o",
   required="True",
   help="Output text file for oncodriveclust detail"
   )

parser$add_argument(
   "--output_plot", "-p",
      required="True",
         help="Output pdf file for oncodriveclust detail"
	    )

args <- parser$parse_args();

###


laml = read.maf(maf = args$input_maf, removeSilent = T, useAll = T)
aacol = 'HGVSp_Short'
if(!is.null(args$aacol)){
aacol = args$aacol
}

if(is.null(args$gene_blacklist)){
	laml.sig = oncodrive(maf =laml, AACol = aacol, pvalMethod = 'zscore',minMut = args$min_mut)
	write.table(laml.sig,file=args$output_detail, quote=FALSE)
	pdf(args$output_plot)
	plotOncodrive(res=laml.sig,fdrCutOff=0.05,useFraction=TRUE)
	dev.off()
}else{
	
	all_genes <- read.table(args$gene_blacklist, stringsAsFactors=FALSE)[,1]

	laml.sig = oncodrive(maf =laml, AACol = aacol, pvalMethod = 'zscore',minMut = args$min_mut,ignoreGenes=all_genes)
	write.table(laml.sig,file=args$output_detail, quote=FALSE)
	pdf(args$output_plot)
	plotOncodrive(res=laml.sig,fdrCutOff=0.05,useFraction=TRUE)
	dev.off()
}