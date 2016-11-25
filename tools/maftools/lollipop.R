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
   "--aacol", "-ac",
   default=0,
   help="Optionally provide the name of the column that contains the amino acid annotation in your MAF file");

parser$add_argument(
   "--output_directory", "-o",
   required="True",
   help="Output directory to store gene plots"
   )

args <- parser$parse_args();

###

all_genes <- read.table(args$gene_list, stringsAsFactors=FALSE)[,1]
mutmaf = read.maf(maf = args$input_maf, removeSilent = T, useAll = T)
count = 1
for (gene in all_genes) {
    output_name = paste(args$output_directory, paste(gene, "_", count, ".pdf", sep=""), sep="/")
    #clean up some oddities in formatting that disagree with MAFtools
    mut = mutmaf@data
    mut$HGVSp_Short = gsub("\\?","*",mut$HGVSp_Short)
    mut$HGVSp_Short = gsub("fs.+","fs",mut$HGVSp_Short)
    mutmaf@data = mut
    pdf(output_name)
    if(args$aacol !=0){
    lollipopPlot(maf = mutmaf, gene = gene, AACol = args$aacol, labelPos = 'all')
    }
    else{
    res = tryCatch({
    lollipopPlot(maf = mutmaf, gene = gene, AACol = 'HGVSp_Short', labelPos = 'all')
    
	}, warning = function(war){
	print("...")
	}, error = function(err){
	print(paste("skipping",gene,"due to lack of protein annotation"))
	})
    }
    dev.off()
    count <- count + 1
    }
