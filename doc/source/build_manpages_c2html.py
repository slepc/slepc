#!/usr/bin/env python3
""" Build and place the manual pages (API) (as .md files) and html source (as .html files)"""

import os
import errno
import subprocess
import shutil
import argparse
import re

rawhtml = ['include', 'src']
#petsc_arch = 'arch-docs'

def pre(doctext,slepc_dir,srcdir,outdir):
    """ Operations to provide data for SLEPc manual pages and c2html files. """
    import time
    import build_man_pages

    x = time.clock_gettime(time.CLOCK_REALTIME)
    print('============================================')
    print('Building all manual pages')
    build_man_pages.main(slepc_dir,srcdir,doctext)
    print("Time: "+str(time.clock_gettime(time.CLOCK_REALTIME) - x))
    print('============================================')

    import build_man_examples_links
    x = time.clock_gettime(time.CLOCK_REALTIME)
    print('============================================')
    print('Building manual page links to tutorials')
    build_man_examples_links.main(slepc_dir,srcdir)
    print("Time: "+str(time.clock_gettime(time.CLOCK_REALTIME) - x))
    print('============================================')

    import build_man_impls_links
    x = time.clock_gettime(time.CLOCK_REALTIME)
    print('============================================')
    print('Building manual page links to implementations')
    build_man_impls_links.main(slepc_dir,srcdir)
    print("Time: "+str(time.clock_gettime(time.CLOCK_REALTIME) - x))
    print('============================================')

    import build_man_index
    x = time.clock_gettime(time.CLOCK_REALTIME)
    print('============================================')
    print('Building manual page indices')
    build_man_index.main(slepc_dir,srcdir)
    print("Time: "+str(time.clock_gettime(time.CLOCK_REALTIME) - x))
    print('============================================')

def post(c2html,mapnames,slepc_dir,srcdir,outdir):
    """ Operations to provide data for SLEPc manual pages and c2html files. """
    import time

#    if not os.path.isfile(os.path.join(slepc_dir, "configure.log")): raise Exception("Expected SLEPc configuration not found")
#    c2html = shutil.which('c2html')
    if not c2html:
      with open(os.path.join(slepc_dir,petsc_arch,'lib','slepc','conf','slepcvariables')) as f:
        c2html = [line for line in f if line.find('C2HTML ') > -1]
        c2html = re.sub('[ ]*C2HTML[ ]*=[ ]*','',c2html[0]).strip('\n').strip()
    print('Using C2HTML:', c2html)
#    mapnames = shutil.which('mapnames')
    if not mapnames:
      with open(os.path.join(slepc_dir,petsc_arch,'lib','slepc','conf','slepcvariables')) as f:
        mapnames = [line for line in f if line.find('MAPNAMES ') > -1]
        mapnames = re.sub('[ ]*MAPNAMES[ ]*=[ ]*','',mapnames[0]).strip('\n').strip()
    print('Using MAPNAMES:', mapnames)
    import build_c2html
    x = time.clock_gettime(time.CLOCK_REALTIME)
    print('============================================')
    print('Building c2html')
    build_c2html.main(slepc_dir,srcdir,outdir,c2html,mapnames)
    print("Time: "+str(time.clock_gettime(time.CLOCK_REALTIME) - x))
    print('============================================')
