#!/usr/bin/env python3
""" Loops through all the source using doctext to generate the manual pages"""

import os
import re
import subprocess
import pathlib

# directories skipped when walking source trees for manual pages
_SKIP_DIRS = ['tests', 'tutorials', 'doc', 'output', 'ftn-custom', 'ftn-auto', 'ftn-mod', 'binding', 'config', 'lib', '.git', 'share', 'systems']

def findlmansec(file):
    mansec = None
    submansec = None
    with open(file) as mklines:
      #print(file)
      submansecl = [line for line in mklines if (line.find('SUBMANSEC') > -1 and line.find('BFORT') == -1)]
      if submansecl:
        submansec = re.sub(r'[ ]*/\* [ ]*SUBMANSEC[ ]*=[ ]*','',submansecl[0]).strip('\n').strip('*/').strip()
        if submansec == submansecl[0].strip('\n'):
          submansec = re.sub('SUBMANSEC[ ]*=[ ]*','',submansecl[0]).strip('\n').strip()
        #print(':SUBMANSEC:'+submansec)
        return submansec
    with open(file) as mklines:
      mansecl = [line for line in mklines if line.startswith('MANSEC')]
      if mansecl:
        mansec = re.sub('MANSEC[ ]*=[ ]*','',mansecl[0]).strip('\n').strip()
        #print(':MANSEC:'+mansec)
        return mansec
    return None

def processdir_batched(slepc_dir, build_dir, dir, doctext):
  '''Runs doctext on batches of source files in the directory'''
  lmansec = None
  if os.path.isfile(os.path.join(dir,'makefile')):
    lmansec = findlmansec(os.path.join(dir,'makefile'))

  batches = []
  for file in os.listdir(dir):
    llmansec = lmansec
    if os.path.isfile(os.path.join(dir,file)) and pathlib.Path(file).suffix in ['.c', '.cxx', '.h', '.cu', '.cpp', '.hpp']:
      if not llmansec:
        llmansec = findlmansec(os.path.join(dir,file))
        if not llmansec: continue
      if not os.path.isdir(os.path.join(build_dir,'manualpages',llmansec)): os.mkdir(os.path.join(build_dir,'manualpages',llmansec))
      if batches and batches[-1][0] == llmansec:
        batches[-1][1].append(file)
      else:
        batches.append((llmansec,[file]))

  numberErrors = 0
  for llmansec,files in batches:
    command = [doctext,
               '-myst',
               '-mpath',    os.path.join(build_dir,'manualpages',llmansec),
               '-heading',  'SLEPc',
               '-defn',     os.path.join(build_dir,'manualpages','doctext','myst.def'),
               '-indexdir', '../'+llmansec,
               '-index',    os.path.join(build_dir,'manualpages','manualpages.cit'),
               '-locdir',   dir[len(slepc_dir)+1:]+'/',
               '-Wargdesc', os.path.join(build_dir,'manualpages','doctext','doctextcommon.txt')]
    sp = subprocess.run(command + files, cwd=dir, capture_output=True, encoding='UTF-8', check=True)
    if sp.stdout and sp.stdout.find('WARNING') > -1:
      print(sp.stdout)
      numberErrors = numberErrors + 1
    if sp.stderr and sp.stderr.find('WARNING') > -1:
      print(sp.stderr)
      numberErrors = numberErrors + 1
  return numberErrors


def main(slepc_dir, srcdir, doctext, use_batch=True):
  # generate the .md files for the manual pages from all the SLEPc source code
  try:
    os.unlink(os.path.join(srcdir,'manualpages','manualpages.cit'))
  except:
    pass
  numberErrors = 0
  skip_dirs = _SKIP_DIRS + [os.environ.get('PETSC_ARCH', 'arch-docs')]
  for dirpath, dirnames, filenames in os.walk(os.path.join(slepc_dir),topdown=True):
    dirnames[:] = [d for d in dirnames if d not in skip_dirs and not d.startswith('arch')]
    numberErrors = numberErrors + processdir_batched(slepc_dir,srcdir,dirpath,doctext)
  if numberErrors:
    raise RuntimeError('Stopping document build since errors were detected in generating manual pages')

  # generate list of all manual pages
  with open(os.path.join(srcdir,'manualpages','htmlmap'),mode='w') as map:
    with open(os.path.join(srcdir,'manualpages','manualpages.cit')) as cit:
      map.write(re.sub(r'man\+../','man+manualpages/',cit.read()))
    with open(os.path.join(srcdir,'manualpages','mpi.www.index')) as mpi:
      map.write(mpi.read())
