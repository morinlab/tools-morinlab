# Cancer Genomics Toolkit for Galaxy
Tool Shed repositories maintained and developed by the Morin Lab

# Overview
This repository hosts all the Galaxy tools, dependency packages and workflows created and/or extended by the Morin laboratory for the Galaxy Cancer Toolkit. 

# Installation
There are several options for installing the tools. The two options we foresee being the likely use cases are detailed below. 

## Option 1: Using the Galaxy Test Toolshed and Your Local Galaxy Installation
These have *largely* been deposited into the Galaxy test toolshed. We are in the process of finalizing this step and appreciate hearing about any lingering issues. One known issue with some of the tools is that they may point to dependencies in another toolshed for legacy reasons. If these are found please let us know and we will correct them. Also, some of the tools do not automatically install dependencies for various reasons (some due to licensing/availability of the software they wrap, others due to issues installing the dependencies on certain platforms).

## Option 2: Use the Docker image
We **strongly** recommend you try using our Dockerfile if you have the capability of running Docker. This will allow you to reproduce a Docker image with the tools described in the manuscript already installed. In principle, this should be a relatively painless way for you to get your hands on the tools without having to set up your own Galaxy server and get the dependencies right as that process can be difficult if not impossible on certain systems. 
