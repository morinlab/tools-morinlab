#!/usr/bin/R

# Script to Run Sequenza Pipeline in Galaxy
library(sequenza);

# INPUT FILE
args=(commandArgs(TRUE));
input.file <- args[1];

ploidy <- args[4];

cellularity <-args[5];

# STEP ONE
extract.data <- sequenza.extract(
file=input.file,
gz=TRUE
);

# STEP TWO
fit.data <- sequenza.fit(
extract.data
);

# STEP THREE
if(length(args)>3){

results.data <- sequenza.results(extract.data, cellularity=cellularity, ploidy=ploidy,out.dir = args[3],sample.id = args[2]);

} else{
  results.data <- sequenza.results(
extract.data,
  fit.data,
out.dir = args[3],
  sample.id = args[2]
);
}
