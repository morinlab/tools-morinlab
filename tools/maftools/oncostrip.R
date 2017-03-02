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
   "--output_plot", "-o",
   required="True",
   help="Output directory to store gene plots"
   )

parser$add_argument(
   "--annotation_file", "-a",
         help="Annotations for samples"
	    )

parser$add_argument(
   "--sort_by", "-s",
   required="True",  
         help="How to sort the matrix before visualizing, allowed options are 'anno' and 'gene'"
	)

args <- parser$parse_args();

###

sort = FALSE
sort_by_anno = FALSE
if(args$sort_by == "anno"){
sort_by_anno=TRUE
print(paste("sort_by_anno",sort_by_anno))
}else if(args$sort_by == "gene"){
sort=TRUE
}else{


}

###

    

all_genes <- read.table(args$gene_list, stringsAsFactors=FALSE)[,1]
laml = read.maf(maf = args$input_maf, removeSilent = T, useAll = T)
count = 1
pdf(args$output_plot,width=20,height=20)
    
        genes.have = laml@gene.summary[,Hugo_Symbol]
	    all_genes = all_genes[all_genes %in% genes.have]


        #oncostrip(maf = laml, gene = all_genes, showTumorSampleBarcodes= FALSE, sort=sort,sort_by_anno=sort_by_anno)
	    if(is.null(args$annotation_file)){
	        oncostrip(maf = laml, gene = all_genes, showTumorSampleBarcodes= FALSE, sort=sort,sort_by_anno=sort_by_anno)
		    }else{
		        df = as.data.frame(read.table(args$annotation_file,sep="\t"))
			        oncostrip(maf = laml, gene = all_genes, showTumorSampleBarcodes= FALSE, sort=sort,sort_by_anno=sort_by_anno,annotation=df)

    }
 dev.off()
