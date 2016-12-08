#!/usr/bin/env python
import sys
import os
from bioblend import galaxy
admin_email = os.environ.get('GALAXY_DEFAULT_ADMIN_USER', 'admin@galaxy.org')
admin_pass = os.environ.get('GALAXY_DEFAULT_ADMIN_PASSWORD', 'admin')
url = "http://localhost:8080"
skip = True

gi = galaxy.GalaxyInstance(url=url, email=admin_email, password=admin_pass)
wf = galaxy.workflows.WorkflowClient(gi)
for arg in sys.argv:
    if skip:
        skip = False
        continue
    infile = arg
    print "importing %s" % infile
    wf.import_workflow_from_local_path(infile)
