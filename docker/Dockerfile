FROM bgruening/galaxy-stable
MAINTAINER Ryan Morin, rdmorin@sfu.ca

RUN apt-get update && apt-get install rsync && apt-get install bc && apt-get install libxp6 libxmu6 #a dependency needed by Strelka installer and GISTIC

RUN add-tool-shed --url 'https://testtoolshed.g2.bx.psu.edu/' --name 'Test Tool Shed'

ENV GALAXY_CONFIG_BRAND Cancer Genomics Galaxy Flavour
ENV GALAXY_CONFIG_ADMIN_USERS  rdmorin@sfu.ca,admin@galaxy.org
ENV GALAXY_CONFIG_LIBRARY_IMPORT_DIR /imported_data/

# Expose port 80 (webserver), 21 (FTP server), 8800 (Proxy)
EXPOSE :80
EXPOSE :21
EXPOSE :8800

RUN sudo apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E084DAB9
RUN echo "deb https://cran.rstudio.com/bin/linux/ubuntu trusty/" >> /etc/apt/sources.list
RUN apt-get update
RUN sh -c '/bin/echo -e "n\n" | apt-get --assume-yes install r-base-dev' #workaround for interactive installer 

RUN apt-get --assume-yes install libcurl4-gnutls-dev #needed for maftools dependencies to install correctly, namely RCurl
RUN apt-get --assume-yes install libgd-dev
RUN apt-get --assume-yes install libfile-copy-recursive-perl
RUN apt-get --assume-yes install Xvfb #needed by IGV screenshot tool
RUN apt-get --assume-yes install libxml2-dev #another dependency for maftools
RUN apt-get --assume-yes install pypy #needed for current implementation of Sequenza to run 
RUN Rscript -e 'install.packages("RCurl",repos="https://cloud.r-project.org/")'
#MAFtools dependency installation is simpler if done outside of Galaxy
RUN Rscript -e 'install.packages("XML",repos="https://cloud.r-project.org/")'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R"); biocLite("ComplexHeatmap"); biocLite("VariantAnnotation"); biocLite("Biostrings")'
RUN Rscript -e 'install.packages("devtools",repos="https://cloud.r-project.org/")'
RUN Rscript -e 'library("devtools"); install_github(repo = "PoisonAlien/maftools");'
RUN Rscript -e 'install.packages("argparse",repos="https://cloud.r-project.org/")'


# Adding the tool definitions to the container
ADD cancer_tools.yml $GALAXY_ROOT/cancer_tools.yml

#add the files that will later be loaded into a history for testing (optional)
ADD https://github.com/morinlab/tools-morinlab/raw/master/test-data/strelka_test_data/test.tumour.bam $GALAXY_ROOT/test.tumour.bam
ADD https://github.com/morinlab/tools-morinlab/raw/master/test-data/strelka_test_data/test.tumour.bam $GALAXY_ROOT/test.normal.bam
ADD https://github.com/morinlab/tools-morinlab/raw/master/test-data/strelka_test_data/test.fa $GALAXY_ROOT/test.fa

# Install the tools
RUN install-tools $GALAXY_ROOT/cancer_tools.yml

#retrieve the custom workflows from github and install them
ADD https://github.com/morinlab/tools-morinlab/raw/master/workflows/Galaxy-Workflow-cnv_merge.ga $GALAXY_ROOT/Galaxy-Workflow-cnv_merge.ga
ADD https://github.com/morinlab/tools-morinlab/raw/master/workflows/Galaxy-Workflow-SignificantGeneDiscovery.ga $GALAXY_ROOT/Galaxy-Workflow-SignificantGeneDiscovery.ga
ADD https://github.com/morinlab/tools-morinlab/raw/master/workflows/Galaxy-Workflow-IGVScreenshotVariants.ga $GALAXY_ROOT/Galaxy-Workflow-IGVScreenshotVariants.ga
ADD https://github.com/morinlab/tools-morinlab/raw/master/workflows/Galaxy-Workflow-MutClusterDiscovery.ga $GALAXY_ROOT/Galaxy-Workflow-MutClusterDiscovery.ga
ADD https://github.com/morinlab/tools-morinlab/raw/master/workflows/Galaxy-Workflow-Ensemble_caller_history.ga $GALAXY_ROOT/Galaxy-Workflow-Ensemble_caller_history.ga
ADD https://github.com/morinlab/tools-morinlab/raw/master/docker/import_workflows.py $GALAXY_ROOT/import_workflows.py

#Warning: If your build fails at this step, try commenting the next two lines out (RUN and ADD)
RUN startup_lite &&  sleep 10 && . $GALAXY_VIRTUAL_ENV/bin/activate &&  python import_workflows.py $GALAXY_ROOT/*.ga

#add script that will populate a history with imported test data (or any other data set desired) after Galaxy is fully launched. 
ADD ./create_and_upload_history.py $GALAXY_ROOT/create_and_upload_history.py

# Mark folders as imported from the host.
VOLUME ["/export/", "/data/", "/var/lib/docker"]
