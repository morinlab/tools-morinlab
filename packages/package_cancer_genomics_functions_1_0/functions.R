base_pairs <- list(
    "A" = "T",
    "C" = "G",
    "G" = "C",
    "T" = "A"
    )

read_biomart_gene_file <- function(biomart_gene_file) {

    data <- read.table(
        file=biomart_gene_file,
        header=F,
       	stringsAsFactors=F
    	)

    data <- data[,c(5,4,3,6)]

    colnames(data) <- c(
        "chr",
        "start",
        "end",
        "gene"
    	)

    data <- data[order(data$end),]
    data <- data[order(data$start),]
    data <- data[order(data$chr),]
    
    return(data)

    }

read_titan_cnv_file <- function(titan_cnv_file, only_major=TRUE, all=FALSE, keep_na=FALSE) {

    data <- read.table(
        file=titan_cnv_file,
        header=T,
        stringsAsFactors=F
    	)

    if ( only_major ) {

        truth <- data[,13] == 1;
        
        if ( keep_na ) {
            truth[is.na(truth)] <- TRUE
            }
        
        else {
            truth[is.na(truth)] <- FALSE
            }
        
        data <- data[truth, ]
        
        }

    if ( !all ) {
        data <- data[,c(2,3,4,1,7,9,10)]

        colnames(data) <- c(
            "chr",
            "start",
            "end",
            "sample",
            "logR",
            "state",
            "copy"
        	)
        }

    return(data)

    }

generate_cohort_cnv_file <- function(list_of_files, method="titan", ...) {

    if ( method == "titan" ) {
        data <- lapply(list_of_files, read_titan_cnv_file, ...);
        data <- Reduce(f=rbind, x=data)
        return(data);
        }

    }

remap_states <- function(cnv_file) {

    copy <- cnv_file$copy;

    new_state <- rep("", length(copy))
    new_state[copy == 0] <- "HOMD";
    new_state[copy == 1] <- "HETD";
    new_state[copy == 2] <- "NEUT";
    new_state[copy == 3] <- "GAIN";
    new_state[copy >= 4] <- "AMP";

    cnv_file$state <- new_state;

    return(cnv_file)

    }

generate_common_genome_frame <- function() {

    data <- data.frame(
        chr = c(1:22, "X", "Y"),
        end = c(
			249250621,
			243199373,
			198022430,
			191154276,
			180915260,
			171115067,
			159138663,
			146364022,
			141213431,
			135534747,
			135006516,
			133851895,
			115169878,
			107349540,
			102531392,
			90354753,
			81195210,
			78077248,
			59128983,
			63025520,
			48129895,
			51304566,
			155270560,
			59373566
			)
        )

    return(data)

    }

calculate_percent_genome_alteration <- function(cnv_frame, genome_frame) {
    
    truth <- !cnv_frame$copy == 2;
    numer <- sum( as.numeric( cnv_frame[truth,]$end - cnv_frame[truth,]$start) );
    denom <- sum( as.numeric( genome_frame$end - 1 ) );
    data <- data.frame(
        sample = cnv_frame$sample[1],
    	pga = numer / denom,
    	stringsAsFactors=F
    	)
    return(data) 

    }

calculate_gene_snv_counts <- function(snv_frame, best=FALSE) {

    data <- snv_frame[,c('gene','patient')];
    binned_data <- table(unique(data)$gene, exclude="Unknown");
    order <- order(binned_data, decreasing=TRUE);
    output <- binned_data[order];
    if (best) {
        output <- output[ output >= mean(output) + 4*sd(output) ];
        }

    output <- data.frame(
        gene = names(output),
        counts = as.numeric(output),
        stringsAsFactors = F
        )

    rownames(output) <- output$gene;

    return(output);    

    }

calculate_patient_snv_counts <- function(snv_frame) {

    data <- snv_frame[,c('patient')];
    binned_data <- table(data);
    order <- order(binned_data, decreasing=TRUE);
    output <- binned_data[order];

    output2 <- data.frame(
        patient = names(output),
        snv_counts = as.numeric(output),
        stringsAsFactors = F
        )

    rownames(output2) <- names(output);

    return(output2);

    }

calculate_patient_snv_frequency <- function(snv_frame) {

    data <- snv_frame[,c('patient','ref','alt','variant_type')];
    data <- data[data$variant_type == "SNP",];

    truth <- data$ref == "A";
    data[truth,]$ref = "T";
    data[truth,]$alt = as.character(base_pairs[data[truth,]$alt]);

    truth <- data$ref == "G";
    data[truth,]$ref = "C";
    data[truth,]$alt = as.character(base_pairs[data[truth,]$alt]);

    data <- data.frame(
        patient = data$patient,
        snv = paste(data$ref, data$alt, sep=">"),
        stringsAsFactors=F
        )

    output <- data.frame(
        patient = rep(unique(data$patient), each=6),
        category = rep(c('C>A', 'C>G', 'C>T', 'T>A', 'T>C', 'T>G'), by=6),
        counts = rep(0, length(unique(data$patient)) * 6),
        stringsAsFactors = F
        );

    tmp <- table(data);

    num_vars <- apply(tmp, 1, sum);

    tmp <- tmp / num_vars;

    for ( patient in rownames(tmp) ) {
        for ( category in colnames(tmp) ) {
            output[output$patient == patient & output$category == category,]$counts <- tmp[patient, category];
            }
        }

    return(output);

    }

read_vcf2maf_maf <- function(vcf2maf_file) {

    data <- read.table(
        vcf2maf_file,
        stringsAsFactors=F,
        header=T,
        comment.char="#",
        sep="\t",
        fill=T,
        quote='"'
        )

    print(names(data))
    output <- data.frame(
        gene = data$Hugo_Symbol,
        patient = data$Tumor_Sample_Barcode,
        chr = data$Chromosome,
        start = data$Start_Position,
        end = data$End_Position,
        variant_class = data$Variant_Classification,
        variant_type = data$Variant_Type,
        ensembl_variant_class = data$Consequence,
        ref = data$Reference_Allele,
        alt = data$Tumor_Seq_Allele2,
        impact = rep("", nrow(data)),
        stringsAsFactors = F
        )

    if (any("IMPACT" == colnames(data))) {
        output$impact = data$IMPACT;
        }

    return(output)

    }

read_gene_list <- function(gene_list_file) {

    data <- read.table(
        gene_list_file,
        stringsAsFactors=F
        )

    output <- data.frame(
        gene = data$V1,
        stringsAsFactors=F
        )

    return(output)

    }


read_mutsigCV_genes <- function(mutsigCV_file) {

    data <- read.table(
        mutsigCV_file,
        stringsAsFactors=F,
        header=T
        )

    data <- data[,c(1,14,15)]

    colnames(data) <- c(
        "gene",
        "p.value",
        "q.value"
        )

    rownames(data) <- data$gene;

    return(data)

    }

read_oncodriveFM_genes <- function(oncodriveFM_file) {

    data <- read.table(
    	oncodriveFM_file,
    	stringsAsFactors=F,
        header=T,
        sep="\t"
    	)

    data <- data[,c(1,2,3)]

    data <- data.frame(
        "gene"=data[,1],
        "p.value"=data[,2],
        "q.value"=data[,3],
        stringsAsFactors=F
    	)

    print(data)

    rownames(data) <- data$gene;

    data <- data[order(data$p.value),]

    return(data)

    }

read_patient_covariate_data <- function(patient_file) {

    data <- read.table(
        patient_file,
        header=T,
        stringsAsFactors=F,
        comment.char=''
        );

    rownames(data) <- data[,1];
    colnames(data)[1] <- "patient";

    return (data);

    }

filter_gene_list <- function(gene_list, q.value=0.05, p.value=1) {
    
    data <- gene_list;
    data <- data[data$q.value <= q.value,];
    data <- data[data$p.value <= p.value,];
    return(data)

    }

join_gene_list <- function(mutsig=NULL, oncodrive=NULL, other=NULL) {

	data <- NULL

	if (!is.null(mutsig)) {
        data <- mutsig;
	    }

	if (!is.null(oncodrive)) {
        if (is.null(data)) {
            data <- oncodrive;
            }
        else {
            data <- rbind(data, oncodrive);
            }
	    }

	if (!is.null(other)) {
        if (is.null(data)) {
            data <- other;
            }
        else {
            data <- rbind(data, other);
            }
	    }

	data <- data[!duplicated(data$gene),]

	return(data)

    }

read_titan_params_file <- function(params_file) {

    data <- read.table(
        params_file,
        fill = TRUE,
        stringsAsFactors=FALSE
        )
   
    return(data)

    }

fetch_validity_index <- function(params_frame) {

    truth <- params_frame[,1] == "S_Dbw" & params_frame[,2] == "validity" &  params_frame[,3] == "index" & params_frame[,4] == "(Both):";

    return(as.numeric(params_frame[truth,5]))
    }

write_cnv_frame <- function(cnv_frame, output, for_circos=FALSE) {

    if ( for_circos ) {
        
        data <- data.frame(
            sample = cnv_frame$sample,
            chr = cnv_frame$chr,
            start = cnv_frame$start,
            end = cnv_frame$end,
            n.mark = rep(1, nrow(cnv_frame)),
            logR = cnv_frame$logR,
            state = cnv_frame$state,
            stringsAsFactors = FALSE
            );

        write.table(
            x = data,
            file = output,
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE,
            sep ="\t"
            );

        }
    else {

        write.table(
            x = cnv_frame,
            file = output,
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t"
            )
        }

    }
