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
    "--gene_mask_list", "-sl",nargs='?',
         help="Input gene blacklist separated by newline (optional)", default=0
);

parser$add_argument(
    "--num_genes", "-n",
        required="True",
	     type="integer",
	         help="Number of genes to include (top N)"
		     );
				      
parser$add_argument(
   "--output_pdf", "-o",
   required="True",
   help="Output directory to store gene plots"
   );

args <- parser$parse_args();

###

if(args$gene_mask_list !=0){
all_genes <- read.table(args$gene_mask_list, stringsAsFactors=FALSE)[,1]
#replace the maf file with a filtered version before loading it
temp.maf = read.csv(args$input_maf,sep="\t",stringsAsFactors = FALSE)
temp.out = temp.maf[!temp.maf[,"Hugo_Symbol"] %in% all_genes,]
write.table(temp.out,file="out.maf",sep="\t",quote=FALSE,row.names=FALSE)
laml = read.maf(maf = "out.maf", removeSilent = T, useAll = T)
pdf(args$output_pdf)
plotmafSummary(maf = laml, rmOutlier = T, addStat = 'median', dashboard = TRUE,top=args$num_genes)
dev.off()
}else{
laml = read.maf(maf = args$input_maf, removeSilent = T, useAll = T)
pdf(args$output_pdf)
plotmafSummary(maf = laml, rmOutlier = T, addStat = 'median', dashboard = TRUE,top=args$num_genes)
dev.off()
}