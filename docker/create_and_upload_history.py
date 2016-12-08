#!/usr/bin/env python
"""
Use the bioblend API to create a fresh history and add a set of files to the history that were imported into the container during the build
Usage: create_and_upload_history.py history_name url1 url2 url3 ...
"""
import sys
from bioblend.galaxy import GalaxyInstance
from bioblend.galaxy.histories import HistoryClient
from bioblend.galaxy.tools import ToolClient

gi = GalaxyInstance(url='http://localhost:80', key='admin')


tc = ToolClient(gi)
lc = HistoryClient(gi)
details = lc.create_history(sys.argv[1])

print "HIST ID: %s" % details["id"]
i = 0
for url in sys.argv:
    url_parts = url.split("/")
    fname = url_parts[-1]
    if i < 2:
        i+=1
        continue
    i+=1
    print "submitting %s as %s" % (url,fname)
    tc.put_url(url,details["id"],file_name=fname)


