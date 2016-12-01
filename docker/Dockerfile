FROM maft_test
MAINTAINER Ryan Morin, rdmorin@sfu.ca

#RUN apt-get update && apt-get install rsync #needed by Strelka installer

#RUN add-tool-shed --url 'http://galaxy-morinlab.mbb.sfu.ca:9009/' --name 'Cancer Genomics Toolshed'
#RUN add-tool-shed --url 'https://testtoolshed.g2.bx.psu.edu/' --name 'Test Tool Shed'



ENV GALAXY_CONFIG_BRAND Cancer Genomics Galaxy Flavour
ENV GALAXY_CONFIG_ADMIN_USERS  rdmorin@sfu.ca,admin@galaxy.org


# Expose port 80 (webserver), 21 (FTP server), 8800 (Proxy)
EXPOSE :80
EXPOSE :21
EXPOSE :8800



# Adding the tool definitions to the container
ADD my_tool_list3.yml $GALAXY_ROOT/my_tool_list3.yml

# Install deepTools
RUN install-tools $GALAXY_ROOT/my_tool_list3.yml
