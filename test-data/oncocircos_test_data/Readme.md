# Oncocircos testing
This tool requires several input files (some may be left empty). Two gene lists allow the user to mask (blacklist.txt) and specifically highlight in color (example_genes_to_label.txt) a subset of genes. An annotation file that provides a unique identifier for each gene alongside genomic coordinates is also required. We have provided one that is easily reproduced using the Ensembl biomart. The main data inputs are a merged set of segmented data similar to that used by GISTIC and the merged MAF-format mutation data from the same cohort. Finally, peak regions deemed significant by GISTIC can also be provided (gistic_sigregions.bed).