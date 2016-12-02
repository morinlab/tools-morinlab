FROM bgruening/galaxy-stable
MAINTAINER Ryan Morin, rdmorin@sfu.ca

RUN apt-get update && apt-get install rsync && apt-get install bc && apt-get install libxp6 libxmu6 #a dependency needed by Strelka installer and GISTIC

RUN add-tool-shed --url 'https://testtoolshed.g2.bx.psu.edu/' --name 'Test Tool Shed'

ENV GALAXY_CONFIG_BRAND Cancer Genomics Galaxy Flavour
ENV GALAXY_CONFIG_ADMIN_USERS  rdmorin@sfu.ca,admin@galaxy.org

# Expose port 80 (webserver), 21 (FTP server), 8800 (Proxy)
EXPOSE :80
EXPOSE :21
EXPOSE :8800

RUN sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
RUN echo "deb https://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list
RUN apt-get update
RUN sh -c '/bin/echo -e "n\n" | apt-get --assume-yes install r-base-dev' #workaround for interactive installer 

RUN apt-get --assume-yes install libcurl4-gnutls-dev #needed for maftools dependencies to install correctly, namely RCurl
RUN apt-get --assume-yes install libxml2-dev #another dependency for maftools
RUN Rscript -e 'install.packages("RCurl",repos="https://cloud.r-project.org/")'
RUN Rscript -e 'install.packages("XML",repos="https://cloud.r-project.org/")'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("ComplexHeatmap"); biocLite("VariantAnnotation"); biocLite("Biostrings")'
RUN Rscript -e 'install.packages("devtools",repos="https://cloud.r-project.org/")'
RUN Rscript -e 'library("devtools"); install_github(repo = "PoisonAlien/maftools");'
RUN Rscript -e 'install.packages("argparse",repos="https://cloud.r-project.org/")'


# Adding the tool definitions to the container
ADD cancer_tools.yml $GALAXY_ROOT/cancer_tools.yml

# Install the tools
RUN install-tools $GALAXY_ROOT/cancer_tools.yml

#set up Matlab dependency and copy required files to Gistic tool directory


RUN cp -R /opt/MATLAB/MATLAB_Compiler_Runtime/ /export/tool_deps/gistic/2.0.22/morinlab/package_gistic_2_0_22/e6bd63a6507a/
