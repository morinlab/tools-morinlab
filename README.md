# Cancer Genomics Toolkit for Galaxy
Tool Shed repositories maintained and developed by the Morin Lab

# Overview
This repository hosts all the Galaxy tools, dependency packages and workflows created and/or extended by the Morin laboratory for the Galaxy Cancer Toolkit. 

# Installation
There are several options for installing the tools. The two options we foresee being the likely use cases are detailed below. 

## Option 1: Using the Galaxy Test Toolshed and Your Local Galaxy Installation
These have *largely* been deposited into the Galaxy test toolshed. We are in the process of finalizing this step and appreciate hearing about any lingering issues. One known issue with some of the tools is that they may point to dependencies in another toolshed for legacy reasons. If these are found please let us know and we will correct them. Also, some of the tools do not automatically install dependencies for various reasons (some due to licensing/availability of the software they wrap, others due to issues installing the dependencies on certain platforms).

## Option 2: Use the Dockerfile provided to build your own image
We **strongly** recommend you try using our [Dockerfile](docker/Readme.md) if you have the capability of running Docker. This will allow you to reproduce a Docker image with the tools described in the manuscript already installed. In principle, this should be a relatively painless way for you to get your hands on the tools without having to set up your own Galaxy server and get the dependencies right as that process can be difficult if not impossible on certain systems. 

## Getting started
We are in the process of putting example data that should demonstrate how to use some of the main tools and workflows provided with this toolkit. Please refer to the documentation in the test-data directory for more details. 

## Known limitations
Please be aware that many of the tools included here have only been tested on a small number of data sets. 

### Not all MAFs are quivalent
We fully expect incompatability with MAF files exported from other pipelines. Ideally, if you are running downstream analyses on your own mutation calls, you should reannotate these. We have packaged the maf2maf tool for this purpose. If you encounter a MAF file that does not work with maf2maf it may be difficult to debug. You can first try paring it down to the minimum required columns for a MAF file. 

###chr1 or 1
We also note that the "chr" prefix is needed/expected by some tools (e.g. GISTIC) but most are either ambivalent or expect no prefix. We suggest that, if possible, you supply reference files and mutation files that lack the chr prefix. Tools that need this are designed to add the prefix if it is lacking. 

###MAF vs MAF
Galaxy has a native MAF file type that represents mutliple sequence alignments. We have not implemented an alternative MAF format yet for the mutation annotation format developed for TCGA. When you upload/create MAF files for this toolkit they will often need to be manually converted to the MAF type before certain tools will run. We expect this minor issue to be fixed soon. 
