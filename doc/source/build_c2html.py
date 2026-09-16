#!/usr/bin/env python3
""" Builds the html files for all the source"""

import os
import re
import subprocess
import pathlib
import concurrent.futures
import multiprocessing
from itertools import chain
from myst_parser.parsers.docutils_ import to_html5_demo

C2HTML_BATCH_SIZE = 32

if __package__:
  from . import build_c2html_file
else:
  import build_c2html_file

def compute_make_np(i):
  '''Number of worker processes to run c2html on'''
  f16 = .80
  f32 = .65
  f64 = .50
  f99 = .30
  if (i<=2):    return 2
  elif (i<=4):  return i
  elif (i<=16): return int(4+(i-4)*f16)
  elif (i<=32): return int(4+12*f16+(i-16)*f32)
  elif (i<=64): return int(4+12*f16+16*f32+(i-32)*f64)
  else:         return int(4+12*f16+16*f32+32*f64+(i-64)*f99)

def main(slepc_dir,srcdir,loc,c2html,mapnames):
  os.chdir(slepc_dir)

  # reformat file that maps manual pages to directory location and add MPI manual pages
  with open('htmlmap.tmp', "w") as fdw, open(os.path.join(srcdir,'manualpages','htmlmap'), "r") as fd:
    fdw.write(fd.read().replace('man+manualpages/','man+HTML_ROOT/manualpages/'))
  with open('htmlmap.tmp', "a") as fdw, open(os.path.join(srcdir,'manualpages','mpi.www.index'), "r") as fd:
    fdw.write(fd.read())

  # walk directories generating list of all source code that needs processing and creating index.html for each directory
  SKIPDIRS = set('public html benchmarks output doc binding config lib bin systems share mpiuni kernels valgrind interfaces data linter'.split())
  SKIPDIRSPREFIX = set('arch- venv- .git'.split())
  SUFFIXES = set('.F90 .F .c .cxx .cpp .h .cu .hpp'.split())
  SUFFIXES_C = set('.c .cxx .cpp .h .cu .hpp'.split())
  SUFFIXES_F = set('.F90 .F'.split())
  sourcefiles = []
  for root, dirs, files in chain.from_iterable(os.walk(path) for path in [slepc_dir]):
    dirs[:] = [d for d in dirs if d not in SKIPDIRS and not any([s for s in SKIPDIRSPREFIX if d.startswith(s)])]
    root = root[len(slepc_dir)+1:]
    if not root: continue
    if not os.path.isdir(os.path.join(loc,root)): os.makedirs(os.path.join(loc,root))
    sourcefiles.extend([os.path.join(root,f) for f in files if any([s for s in SUFFIXES if f.endswith(s)])])

    # create index.html file for each directory
    with open(os.path.join(loc,root,'index.html'),'w') as fdw:
      if root.startswith('src'):

        # get MANSEC from the makefile and copy the MANSEC basic information into the index
        if os.path.isfile(os.path.join(root,'makefile')):
          with open(os.path.join(root,'makefile')) as mklines:
            mansecl = [line for line in mklines if line.startswith('MANSEC')]
            if mansecl:
              mansec = re.sub('MANSEC[ ]*=[ ]*','',mansecl[0]).strip('\n').strip()
              with open(os.path.join('doc','source','manualpages','MANSECHeaders',mansec)) as fd:
                for line in fd:
                  if (line.find('Related Users Manual part:') > -1 or
                      line.find('>Examples</a>') > -1):
                    continue
                  fdw.write(to_html5_demo(line))

      fdw.write('\n<p>\n')

      # TODO: use HTML lists for the list below

      # list examples
      if root.find('/tests') > -1 or root.find('tutorials') > -1:
        fdw.write('\n<p>\nExamples\n<p>')
        examples = {}
        for f in files:
          if any([s for s in SUFFIXES_C if f.endswith(s)]):
            with open(os.path.join(root,f)) as fd:
              examples[f] = ''
              for line in fd:
                l = line.find('char help[] = ')
                if l > -1:
                  s = line.find('\\n')
                  examples[f] = line[l + 15:s]
                  break
        for f in files:
          if any([s for s in SUFFIXES_F if f.endswith(s)]):
            with open(os.path.join(root,f)) as fd:
              examples[f] = ''
              for line in fd:
                l = line.find('Description:')
                if l > -1:
                  examples[f] = line[l + 13:]
                  break
        # simple natural sorting
        examples = dict(sorted(examples.items(), key=lambda i: [int(s) if
                                                           s.isdigit() else
                                                           s.lower() for s in
                                                           re.split(r'(\d+)',
                                                                    i[0])]))
        for f in examples.keys():
          fdw.write('<a href="' + f + '.html">' + f + ': ' + examples[f] + '</a><br>\n')

      # list source code
      else:
        if any([f for f in files if any([s for s in SUFFIXES if f.endswith(s)])]):
          fdw.write('\n<p>\nSource files\n<p>')
        for f in files:
          if any([s for s in SUFFIXES if f.endswith(s)]):
            fdw.write('<a href=\"' + f + '.html\">' + f + '</a><br>\n')

      # list subdirectories
      if dirs:
        fdw.write('\n<p>\nDirectories\n<p>')
      for d in dirs:
        fdw.write('<a href="' + os.path.join(d,'index.html') + '">' + d + '</a><br>\n')

  git_sha = subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).rstrip()
  # Reuse a fixed set of Python interpreters instead of starting one for every source file.
  with concurrent.futures.ProcessPoolExecutor(max_workers=compute_make_np(multiprocessing.cpu_count())) as executor:
    futures = [
      executor.submit(build_c2html_file.main_batch,slepc_dir,loc,git_sha,c2html,mapnames,sourcefiles[i:i + C2HTML_BATCH_SIZE])
      for i in range(0,len(sourcefiles),C2HTML_BATCH_SIZE)
    ]
    for future in concurrent.futures.as_completed(futures):
      future.result()

  os.unlink('htmlmap.tmp')
