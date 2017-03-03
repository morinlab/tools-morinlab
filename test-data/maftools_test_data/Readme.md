# MAFtools suite of Galaxy visualization tools
These files can be used to reproduce a variety of the images generated with our Galaxy tools based on the [MAFtools R package](https://github.com/PoisonAlien/maftools). Upload the gzipped maf file in this directory to a Galaxy instance with the maftools package installed and ensure that you set the data type to "MAF". You will also want to upload the provided blacklist file or provide your own, to suppress the display of unwanted genes in some visualizations. We realize that this clashes with another distinct data type in Galaxy and hope to remedy this ambiguity with future releases. Suggestions and/or pull requests are welcome!

## Reproducing specific figures from our paper
### GeneCloud in Figure 3, panel B
Run the Gene Cloud Plot tool on this data with the MAF and (optionally) blacklist file as input. 

### MAF summary in Figure 3, panel C
Run the Maf Summary Plot tool on the provided MAF file and the blacklist file. 

### Oncoprintplus in Figure 4, panel C
Run the oncoprintplus tool with the provided MAF and the outputs from oncodriveFM and/or MutSig. An example workflow showing how to prepare your input for oncodriveFM is included in the [provided workflow](https://github.com/morinlab/tools-morinlab/tree/master/workflows#using-oncodrivefm-to-detect-significantly-mutated-genes-and-visualizations). This workflow will also generate the Lollipop Plots for each significant gene, such as those examples shown in Panel C of the same figure.



