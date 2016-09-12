suppressMessages(library(BoutrosLab.plotting.general));
suppressMessages(library(dplyr));

### GET ARGUMENTS ##################################################################################

suppressMessages(library(argparse));
parser <- ArgumentParser(description="Create a Gene by Patient Heatmap BPGs multiplot fucntion");

### VARIANT INPUT ARGUMENTS ########################################################################

parser$add_argument(
    "--input_snv", "-snv",
    required="True",
    help="Input SNV files in MAF format"
	);

parser$add_argument(
    "--input_cnv", "-cnv",
    help="Input CNV files in MAF format"
	);


### GENE INPUT ARGUMENTS ###########################################################################

parser$add_argument(
    "--gene_list", "-gl",
    help=""
	);

parser$add_argument(
    "--gene_oncodrive", "-go",
    help=""
	);

parser$add_argument(
    "--gene_mutsig", "-gc",
    help=""
	)

### GENERAL ARGUMENTS ##############################################################################

parser$add_argument(
    "--center_plot",
    choices=c('impact', 'variant_classification'),
    required=TRUE,
    help="What should be displayed in the Center Heatmap?"
	);

parser$add_argument(
    "--cancer_genomics_functions_path", "-cgfp",
    required='True',
    help="Path to cancer genomics functions in R"
    )

parser$add_argument(
	"--output", "-o",
    required="True",
    help="Output Image File"
	)

### PATIENT ARGUMENTS ##############################################################################

parser$add_argument(
    "--patient_covariate_data",
    help="Input Patient Covariate Matrix"
    );

parser$add_argument(
    "--plot_patient_covariate",
    action="append",
    help="Should We Plot the Input Covariate Data?"
    );

parser$add_argument(
    "--plot_patient_snv_counts",
    action='store_true',
    help="Should We Plot Patient SNV Counts?"
	);

parser$add_argument(
    "--plot_patient_snv_distribution",
    action='store_true',
    help="Should We Plot Patient SNV Distribution?"
	);

parser$add_argument(
	"--patient_order",
    action="append",
    help="How Should We Determine Patient Order?"
	);


### GENE ARGUMENTS #################################################################################

parser$add_argument(
    "--plot_gene_snv_counts",
    action="store_true",
    help="Should We Plot Gene SNV Counts?"
	);

parser$add_argument(
	"--plot_gene_mutsig",
	action="store_true",
	help="Should We Plot Gene Significance?"
	);

parser$add_argument(
    "--plot_gene_oncodrive",
    action="store_true",
    help="How Should We Determine Gene Order?"
	);

parser$add_argument(
    "--gene_order",
    choices=c('mutsig', 'oncodrive', 'snv_counts', 'gene_list'),
    help="How Should We Determine Gene Order?"
	);

####################################################################################################

args <- parser$parse_args();

source(paste(args$cancer_genomics_functions_path, "functions.R", sep="/"));

### ORGANIZE GENE DATA #############################################################################

snv_data <- read_vcf2maf_maf(args$input_snv);

snv_data$patient

snv_counts <- calculate_gene_snv_counts(snv_data);

patient_covariate_data <- calculate_patient_snv_counts(snv_data);

patient_snv_distribution <- calculate_patient_snv_frequency(snv_data);

if (!is.null(args$patient_covariate_data)) {
    tmp_patient_covariate_data <- read_patient_covariate_data(args$patient_covariate_data);

    if (!all(rownames(patient_covariate_data) %in% rownames(tmp_patient_covariate_data))) {
        stop("Some patients identified in the MAF file are not present in patient covariate data file");   
        }
    indexes <- which( rownames(patient_covariate_data) %in% rownames(tmp_patient_covariate_data))
    current_names <- colnames(patient_covariate_data);
    patient_covariate_data <- cbind(patient_covariate_data, tmp_patient_covariate_data[rownames(patient_covariate_data)[indexes],-1]);
    colnames(patient_covariate_data) <- c(current_names, colnames(tmp_patient_covariate_data)[-1]);
    }


##################################

gene_list <- NULL;
mutsig <- NULL;
mutsig_sig <- NULL;
oncodrive <- NULL;
oncodrive_sig <- NULL;

if (!is.null(args$gene_list)) {
    gene_list <- read_gene_list(args$gene_list);
    }
if (!is.null(args$gene_mutsig)) {
    mutsig <- read_mutsigCV_genes(args$gene_mutsig);
    mutsig_sig <- filter_gene_list(mutsig);
    };
if (!is.null(args$gene_oncodrive)) {
    oncodrive <- read_oncodriveFM_genes(args$gene_oncodrive);
    oncodrive_sig <- filter_gene_list(oncodrive);
    };

genes_of_interest <- c(
	gene_list$gene, 
	mutsig_sig$gene, 
	oncodrive_sig$gene
	);

if (!is.null(args$gene_list)) {
    genes_of_interest <- gene_list$gene
    }

genes_of_interest <- unique(genes_of_interest);

if (is.null(genes_of_interest)) {
    genes_of_interest <- calculate_gene_snv_counts(snv_data, best=TRUE)$gene
    }

gene_data <- data.frame(
    gene = genes_of_interest,
    user = genes_of_interest %in% gene_list,
    mutsig = genes_of_interest %in% mutsig_sig$gene,
    oncodrive = genes_of_interest %in% oncodrive_sig$gene,
    stringsAsFactors = F
	);

gene_order <- NULL;

if (is.null(args$gene_order)) {
    gene_order <- gene_data$gene;

} else if (args$gene_order == 'mutsig') {
	add_to_end <- NULL;
	if (!all(gene_data$gene %in% mutsig$gene)) {
        warning("mutsigCV input file is missing some genes, assuming a p value of 1");
        add_to_end <- gene_data$gene[!gene_data$gene %in% mutsig$gene]
	    }
    
    gene_order <- c(mutsig$gene[mutsig$gene %in% gene_data$gene], add_to_end);

} else if (args$gene_order == 'oncodrive') {
	add_to_end <- NULL;
	if (!all(gene_data$gene %in% oncodrive$gene)) {
        warning("oncodriveFM input file is missing some genes, assuming a p value of 1");
        add_to_end <- gene_data$gene[!gene_data$gene %in% oncodrive$gene]
	    }
    
    gene_order <- c(oncodrive$gene[oncodrive$gene %in% gene_data$gene], add_to_end);  
    
} else if (args$gene_order == 'snv_counts') {
    add_to_end <- NULL;
    if ( !all(gene_data$gene %in% snv_counts$gene)) {
        warning("snv input file is missing some genes, assuming snv count of 0");
        add_to_end <- gene_data$gene[!gene_data$gene %in% snv_counts$gene]
        }

    gene_order <- c( snv_counts[snv_counts$gene %in% gene_data$gene,]$gene, add_to_end);
    }
    
patient_order <- NULL;



if ( is.null(args$patient_order) ) {
    patient_order <- unique(snv_data$patient);
    patient_covariate_data <- patient_covariate_data[patient_order,];

} else {
    for (covariate in args$patient_order) {
        print(covariate)
        print(colnames(patient_covariate_data))
        current_order <- order(patient_covariate_data[,covariate]);
        patient_covariate_data <- patient_covariate_data[current_order,];
        } 

    patient_order <- rownames(patient_covariate_data);

    }

### GENE PLOTS #####################################################################################

gene_plots <- list();


if (args$plot_gene_snv_counts) {

    plot_data <- data.frame(
        gene = gene_order,
        order = length(gene_order):1,
        counts = log10(snv_counts[gene_order,]$counts)
    	);

    max <- ceiling(max(plot_data$counts))
    middle <- round(median(plot_data$counts))

    plot <- create.barplot(
        data = plot_data,
        formula = order ~ counts,
        plot.horizontal = TRUE,
        xlimits = c(-0.1, max+0.1),
        xat = c(0, middle, max),
        xaxis.lab = c( expression(10^0), bquote(10^.(middle)), bquote(10^.(max)) ),
        yat = 0
       	);

    gene_plots <- c(gene_plots, list(plot));

    }


if (args$plot_gene_mutsig) {

    plot_data <- data.frame(
        gene = gene_order,
        order = length(gene_order):1,
        p.value = -1*log10(mutsig[gene_order,]$p.value),
        q.value = mutsig[gene_order,]$q.value,
        stringsAsFactors=F
    	);

    if (any( is.infinite(plot_data$p.value) ) ) {
        plot_data$p.value[is.infinite(plot_data$p.value)] <- max( c( plot_data$p.value[!is.infinite(plot_data$p.value)], 10) );
        }

    plot_data$p.value[is.na(plot_data$p.value)] <- 1
    plot_data$q.value[is.na(plot_data$q.value)] <- 1
    plot_data$colour = c("black", "grey80")[(plot_data$q.value <= 0.05) + 1]

    max <- ceiling(max(plot_data$p.value))
    middle <- round(median(plot_data$p.value))

    plot <- create.barplot(
        data = plot_data,
        formula = order ~ p.value,
        col = rev(plot_data$colour),
        plot.horizontal = TRUE,
        xlimits = c(-1, max+1),
        xat = c(0, middle, max),
        xaxis.lab = c( expression(10^-0), bquote(10^-.(middle)), bquote(10^-.(max)) ) ,
        yat = 0,
       	);

    gene_plots <- c(gene_plots, list(plot))

    }


if (args$plot_gene_oncodrive) {

    plot_data <- data.frame(
        gene = gene_order,
        order = length(gene_order):1,
        p.value = -1 * log10(oncodrive[gene_order,]$p.value),
        q.value = oncodrive[gene_order,]$q.value,
        stringsAsFactors=F
    	);

    plot_data$p.value[is.na(plot_data$p.value)] <- 1
    plot_data$q.value[is.na(plot_data$q.value)] <- 1
    plot_data$colour = c("black", "grey80")[(plot_data$q.value <= 0.05) + 1]

    max <- ceiling(max(plot_data$p.value))
    middle <- round(median(plot_data$p.value))

    plot <- create.barplot(
        data = plot_data,
        formula = order ~ p.value,
        col = rev(plot_data$colour),
        plot.horizontal = TRUE,
        xlimits = c(-1, max),
        xat = c(0, middle, max),
        xaxis.lab = c( expression(10^-0), bquote(10^-.(middle)), bquote(10^-.(max)) ),
        yat = 0
       	);

    gene_plots <- c(gene_plots, list(plot))   

    }


### PATIENT PLOTS ##################################################################################

patient_plots <- list();

if (args$plot_patient_snv_counts) {

	plot_data <- data.frame(
        patient = patient_order,
        order = 1:length(patient_order),
        counts = log10(patient_covariate_data[patient_order,]$snv_counts)
        )

    middle <- round(median(plot_data$counts))
    max <- ceiling(max(plot_data$counts))

	plot <- create.barplot(
        data = plot_data,
        formula = counts ~ order,
        ylimits = c(-0.1, max+0.1),
        yaxis.lab = c(expression(10^0), bquote(10^.(middle)), bquote(10^.(max)) ),
        yat = c(0,middle,max),
        xat=0
		);

    patient_plots <- c(patient_plots, list(plot));

    }

if (args$plot_patient_snv_distribution) {

    patient_to_order <- data.frame(
        order = 1:length(patient_order)
    	);

    rownames(patient_to_order) <- patient_order;

    plot_data <- data.frame(
        patient =  patient_snv_distribution$patient,
        category = patient_snv_distribution$category,
        counts = patient_snv_distribution$counts,
        order = patient_to_order[patient_snv_distribution$patient,],
        colours = rep(c("#0e8efa", "#010101","#fe0000", "#bfbfbf", "#00f100", "#ffbfcb"), by=6),
        stringsAsFactors=F
    	);
    
    plot <- create.barplot(
        data = plot_data,
        formula = counts ~ order,
        groups = category,
        stack = TRUE,
        xat = 0,
        yaxis.lab = c("0.0", "0.2", "0.4", "0.6", "0.8", "1.0"),
        ylimits = c(-0.1, 1.1),
        yat = c(0, .2, .4, .6, .8, 1.0),
        col = plot_data$colours
    	);

    patient_plots <- c(patient_plots, list(plot));

    }

### PATIENT COVARIATE PLOT #########################################################################

patient_covariate_plots <- list();
patient_covariate_legend <- list();

if (!is.null(args$plot_patient_covariate)) {

    cov.data <- patient_covariate_data[patient_order,c('snv_counts', args$plot_patient_covariate)]; 

    cov.colors <- c();
    for (cov in args$plot_patient_covariate) {
        tmp <- paste0(cov, "_col");
        if ( any(tmp == colnames(patient_covariate_data)) ) {
            tmp.colors <- unique(patient_covariate_data[,tmp]);
            for (col_index in (length(cov.colors)+1):((length(cov.colors)+1)+length(tmp.colors))) {
                col <- tmp.colors[col_index - length(cov.colors)];
                cov.data[which(col == patient_covariate_data[,tmp]), cov] <- col_index;
                }
            }
        cov.colors <- c(cov.colors, tmp.colors);
        }
    cov.data <- cov.data[,-1];

    covariates <- create.heatmap(
        x = matrix(as.numeric(as.matrix(cov.data)), ncol=length(args$plot_patient_covariate)),
        clustering.method = 'none',
        colour.scheme = as.vector(cov.colors),
        total.colours = length(cov.colors)+1,
        row.colour = 'white',
        col.colour = 'white',
        grid.row = TRUE,
        grid.col = TRUE,
        yaxis.lab = args$plot_patient_covariate,
        yat = 1:length(args$plot_patient_covariate),
        xaxis.lab = patient_order,
        xat = 1:length(patient_order),
        yaxis.tck = 0,
        print.colour.key = FALSE
        );

    patient_covariate_plots <- c(patient_covariate_plots, list(covariates));

    for ( cov in args$plot_patient_covariate) {

        legend <- list(
            legend = list(
                colours = rev(unique(patient_covariate_data[, paste0(cov, "_col")])),
                labels = rev(unique(patient_covariate_data[, cov])),
                border = 'white',
                title = cov,
                pch = 15
                )
            );

        patient_covariate_legend <- c(patient_covariate_legend, legend);

        }

    }

### CENTRAL PLOT ###################################################################################

central_plot <- NULL;

if ( args$center_plot == 'impact' ) {

    plot_data <- matrix(
        data = rep(0, length(patient_order) * length(gene_order)),
        nrow = length(gene_order)
    	)

    rownames(plot_data) <- rev(gene_order);
    colnames(plot_data) <- patient_order;

    for (patient in patient_order) {
        for (gene in gene_order) {
        	truth <- snv_data$patient == patient & snv_data$gene == gene;
            if (any(truth)) {
                if (any(snv_data$impact[truth] == "HIGH")) {
                    plot_data[gene, patient] <- 3;
                } else if (any(snv_data$impact[truth] == "MODERATE")) {
                	plot_data[gene, patient] <- 2;
                } else if (any(snv_data$impact[truth] == "LOW")) {
                	plot_data[gene, patient] <- 1;
                    }
                }
            }
        }

    gene.labs <- rownames(plot_data);
    patient.labs <- colnames(plot_data);
    if (!is.null(args$plot_patient_covariate)) {
        patient.labs <- rep("", length(patient.labs));
        }

    plot <- create.heatmap(
	    t(plot_data),
	    clustering.method = "none",
	    colour.scheme = c('grey95',c('#ffeda0', '#feb24c','#f03b20')),
	    total.colours = 5,
	    xaxis.tck = 0,
	    yaxis.tck = 0,
	    print.colour.key = F,
        xaxis.lab = patient.labs,
        yaxis.lab = gene.labs,
	    row.colour = 'white',
	    col.colour = 'white',
	    grid.row = TRUE,
	    grid.col = TRUE
	    );

    central_plot <- list(plot);

    }

which_present <- rep(FALSE, 9)

if (args$center_plot == 'variant_classification') {

    plot_data <- matrix(
        data = rep(0, length(patient_order) * length(gene_order)),
        nrow = length(gene_order)
    	)

    rownames(plot_data) <- rev(gene_order);
    colnames(plot_data) <- patient_order;

	order.of.variants <- data.frame(
	    Truncation = 'Truncation',
	    Splicing = 'Splicing',
	    Missense = 'Missense',
	    UTR = 'UTR',
        Non_coding = "Non_coding",
        Regulatory = "Regulatory",
        miRNA = "miRNA",
        Synonymous = "Synonymous",
	    stringsAsFactors = F
		)

    categories <- recode(
        snv_data[,]$ensembl_variant_class,
        # Truncating variants
        stop_gained = "Truncation",
        frameshift_variant = "Truncation",
        stop_lost = "Truncation",
        start_lost = "Truncation",
        incomplete_terminal_codon_variant = "Truncation",
        transcript_ablation = "Truncation",
        # Missense(-ish) variants
        missense_variant = "Missense",
        inframe_deletion = "Missense",
        inframe_insertion = "Missense",
        protein_altering_variant = "Missense",
        # UTR Variants
        `5_prime_UTR_variant` = "UTR",
        `3_prime_UTR_variant` = "UTR",
        # Splicing variants
        splice_acceptor_variant = "Splicing",
        splice_donor_variant = "Splicing",
        splice_region_variant = "Splicing",
        # Non-coding variants
        non_coding_transcript_exon_variant = "Non_coding",
        intron_variant = "Non_coding",
        upstream_gene_variant = "Non_coding",
        downstream_gene_variant = "Non_coding",
        intergenic_variant = "Non_coding",
        # Regulatory variants
        TF_binding_site_variant = "Regulatory",
        regulatory_region_variant = "Regulatory",
        # Synonymous variants
        synonymous_variant = "Synonymous",
        coding_sequence_variant = "Synonymous",
        stop_retained_variant = "Synonymous",
        # miRNA variants
        mature_miRNA_variant = "miRNA",
        # Default
        .default = NA_character_)

	for( gene in gene_order ) {
	    for( patient in patient_order ) {
	    	truth <- snv_data$patient == patient & snv_data$gene == gene;
            if (gene == "SRSF7" && any(truth)) {
                }
	        if (any(truth)) {
	        	vars <- colnames(order.of.variants) %in% categories[truth]
	        	if(any(vars)) {
	                plot_data[gene, patient] <- which(vars)[1]
	                which_present[which(vars)[1]] <- TRUE
                    }
                else {
                    which_present[9] <- TRUE
                    }
	            }
	        }
	    }

    col_one <- "#512C6F";
    col_two <- "#0F6A99";
    col_thr <- "#46823C";
    col_fou <- "#B367A7";
    col_fiv <- "#64B4D5";
    col_six <- "#F7D72E";
    col_sev <- "#EF922A";
    col_eig <- "#B12025";

    gene.labs <- rownames(plot_data);
    patient.labs <- colnames(plot_data);
    if (!is.null(args$plot_patient_covariate)) {
        patient.labs <- rep("", length(patient.labs));
        }

    plot <- create.heatmap(
	    t(plot_data),
	    clustering.method = "none",
	    colour.scheme = c('grey95', c(col_one, col_two, col_thr, col_fou, col_fiv, col_six, col_sev, col_eig) ),
	    total.colours = 10,
	    xaxis.tck = 0,
	    yaxis.tck = 0,
	    print.colour.key = F,
        xaxis.lab = patient.labs,
        yaxis.lab = gene.labs,
	    row.colour = 'white',
	    col.colour = 'white',
	    grid.row = TRUE,
	    grid.col = TRUE
	    );

    central_plot <- list(plot);

    }

### LEGENDS ########################################################################################

legends <- list();

if (args$plot_patient_snv_distribution ) {

	legend <- list(
		legend = list(
		    colours = rev(c("#0e8efa", "#010101","#fe0000", "#bfbfbf", "#00f100", "#ffbfcb" )),
		    labels = rev(expression(C > A,C > G,C > T, T > A, T > C, T > G)),
		    border = 'white',
		    title = 'Mutation Frequencies',
		    pch = 15
		    )
		);

	legends <- c(legends, legend);

    }

if (args$plot_gene_mutsig || args$plot_gene_oncodrive) {
	
	legend <- list(      
	    legend = list(
	    	colours = c('grey80', 'black'),
	    	labels = expression(Q <= 10^-1, Q > 10^-1),
	    	border = 'white',
	    	title = 'Q-value Groupings',
	    	pch=15
	    	)
	    )

    legends <- c(legends, legend);

    }

if ( args$center_plot == 'impact' ) {

	legend <- list(   
	    legend = list(
	        colours = c('grey95',c('#ffeda0', '#feb24c','#f03b20')),
	        labels = c("None", "Low", "Moderate", "High"),
	        border = 'white',
	        title = 'SNV Impact',
	        pch = 15
	        )
	    )

    legends <- c(legends, legend);

    }

if ( args$center_plot == 'variant_classification' ) {

    col_one <- "#512C6F";
    col_two <- "#0F6A99";
    col_thr <- "#46823C";
    col_fou <- "#B367A7";
    col_fiv <- "#64B4D5";
    col_six <- "#F7D72E";
    col_sev <- "#EF922A";
    col_eig <- "#B12025";

	legend <- list(   
	    legend = list(
	        colours = c(c(col_one, col_two, col_thr, col_fou, col_fiv, col_six, col_sev, col_eig), 'grey95' )[which_present],
	        labels = c("Truncation", "Splicing", "Missense", "UTR", "Non-coding", "Regulatory", "miRNA", "Synonymous", "None")[which_present],
	        border = 'white',
	        title = 'SNV Category',
	        pch = 15
	        )
	    )

    legends <- c(legends, legend);

    }

if (length(patient_covariate_legend) > 0) {

    legends <- c(legends, patient_covariate_legend);

    }

legend_grob <- legend.grob(
    legends = legends,
    title.just = 'left',
    label.cex = 1.0,
    title.cex = 1.0
    );

### MULTIPLOT ######################################################################################

plots <- c(
    patient_covariate_plots,
    central_plot,
    gene_plots,
    patient_plots
	);

plot.layout <- c(1 + length(gene_plots), 1 + length(patient_plots));
layout.skip <- c(F, rep(F, length(gene_plots)), rep(c(F, rep(T, length(gene_plots))), length(patient_plots)));
panel.heights <- c(rep(3, length(patient_plots)),20);
panel.widths <- c(20,rep(3, length(gene_plots)));
if (length(patient_covariate_plots) > 0) {
    plot.layout[2] <- plot.layout[2] + 1;
    layout.skip <- c(F, rep(T, length(gene_plots)), layout.skip);
    panel.heights <- c(panel.heights, 0.5 * length(args$plot_patient_covariate) );
    }

create.multiplot(
    filename = args$output,
    res = 300,
    plot.objects = plots,
    plot.layout = plot.layout,
    layout.skip = layout.skip,
    panel.heights = panel.heights,
    panel.widths = panel.widths,
    height = 16,
    width = 16,
    x.spacing = -.8,
    y.spacing = -.33 * max(nchar(patient_order)),
    left.padding = 4,
    xaxis.cex = .65,
    yaxis.cex = .9,
    xaxis.rot = 90,
    retrieve.plot.labels=T,
    legend = list(
		right = list(
			x = 0.10,
			y = 0.50,
			fun = legend_grob
			)
		),

    print.new.legend= TRUE

	);
