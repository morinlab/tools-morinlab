library(argparse)

parser <- ArgumentParser(description="Select the best TITAN cluster");

parser$add_argument(
    '--input_params', "-ip",
    nargs="+",
    required="True",
    help="Input parameter files from TITAN"
    );

parser$add_argument(
    '--input_segments', "-is",
    nargs="+",
    required="True",
    help="Input segment files from TITAN"
    );

parser$add_argument(
    "--method", "-m",
    choices=c('validity_index', 'validity_index_and_pga'),
    required="True",
    help="Algorithm to use for cluster selection"
    );

parser$add_argument(
    "--percent_genome_alteration", "-pga",
    metavar="PGA",
    type="double",
    default="0.8",
    help="Maximum percent genome alteration to accept in a cluster"
    );

parser$add_argument(
    "--output_param", "-op",
    required="True",
    help="Output file to store optimal parameter file"
    );

parser$add_argument(
    "--output_segment", "-os",
    required="True",
    help="Output file to store optimal segment file"
    );

parser$add_argument(
    "--cancer_genomics_functions_path", "-cgfp",
    required="True",
    help="Path to cancer genomics functions in R"
    )

args <- parser$parse_args();

source(paste(args$cancer_genomics_functions_path, "functions.R", sep="/"));

current_index <- 0;
current_val_index <- 1;
genome_frame <- generate_common_genome_frame();
for (i in 1:length(args$input_params)) {
    if (args$method == "validity_index_and_pga") {
        titan_segs <- read_titan_cnv_file(args$input_segments[i]);
        pga <- calculate_percent_genome_alteration(titan_segs, genome_frame);
        if (pga >= args$percent_genome_alteration) {
            next;
            }
        }

    titan_params <- read_titan_params_file(args$input_params[i]);
    val_index <- fetch_validity_index(titan_params)
    if (current_val_index > val_index) {
        current_val_index <- val_index;
        current_index <- i;
        }
    }

command_params <- paste(
    "cp",
    args$input_params[current_index],
    args$output_param,
    sep=" "
    );

system(command_params);


command_segments <- paste(
    "cp",
    args$input_segments[current_index],
    args$output_segment,
    sep=" "
    );

system(command_params);
