args <- commandArgs(trailingOnly = TRUE)


my_read_table <- function(x) {
    return(read.table(x, stringsAsFactors=F))
    }


data_files <- my_read_table(args[1])
interval_files <- my_read_table(args[2])
byte_files <- my_read_table(args[3])
order <- my_read_table(args[4])[,1]
output <- args[5]


interval_file_list <- apply(interval_files, 1, my_read_table)
byte_file_list <- apply(byte_files, 1, my_read_table)

which_equal <- function(a,b) {
    which(a == b);
    }

find_indexes <- function(element, a_list) {
    a <- lapply(a_list, which_equal, element);
    list_index <- which(!is.na(a==0));
    vector_index <- a[[list_index]];
    return(c(list_index, vector_index));
    }

data <- t(sapply(order, find_indexes, interval_file_list));

print(interval_file_list)
print(byte_file_list)

for ( i in 1:nrow(data)) {

    command <- paste(
        "tail -c +", 
        byte_file_list[[ data[i,1] ]][,1][ data[i,2] ],
        " ",
        data_files[ data[i,1], 1 ],
        " | ", 
        "head -c ", 
        byte_file_list[[ data[i,1] ]][,1][ data[i,2]+1 ] - byte_file_list[[ data[i,1] ]][,1][ data[i,2] ],
        " >> ", 
        output, 
        sep=""
        );

    print(command)
    
    system(command);
    
    }

