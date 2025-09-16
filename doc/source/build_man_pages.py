#!/usr/bin/env python3
""" Loops through all the source using doctext to generate the manual pages"""

import os
import re
import subprocess
import pathlib

def findlmansec(file):
    mansec = None
    submansec = None
    with open(file) as mklines:
#      print('  '+file)
      submansecl = [line for line in mklines if (line.find('SUBMANSEC') > -1 and line.find('BFORT') == -1)]
      if submansecl:
        submansec = re.sub(r'[ ]*/\* [ ]*SUBMANSEC[ ]*=[ ]*','',submansecl[0]).strip('\n').strip('*/').strip()
        if submansec == submansecl[0].strip('\n'):
          submansec = re.sub('SUBMANSEC[ ]*=[ ]*','',submansecl[0]).strip('\n').strip()
#        print('  :SUBMANSEC:'+submansec)
        return submansec
    with open(file) as mklines:
      mansecl = [line for line in mklines if line.startswith('MANSEC')]
      if mansecl:
        mansec = re.sub('MANSEC[ ]*=[ ]*','',mansecl[0]).strip('\n').strip()
#        print('  :MANSEC:'+mansec)
        return mansec
    return None

def processdir(slepc_dir, srcdir, dir, doctext):
  '''Runs doctext on each source file in the directory'''
#  print('Processing dir: '+dir)
  doctext_path = os.path.join(srcdir,'manualpages','doctext')
  lmansec = None
  if os.path.isfile(os.path.join(dir,'makefile')):
    lmansec = findlmansec(os.path.join(dir,'makefile'))

  numberErrors = 0
  for file in os.listdir(dir):
    llmansec = lmansec
    if os.path.isfile(os.path.join(dir,file)) and pathlib.Path(file).suffix in ['.c', '.cxx', '.h', '.cu', '.cpp', '.hpp']:
#      print(' Processing file: '+file)
      if not llmansec:
        llmansec = findlmansec(os.path.join(dir,file))
        if not llmansec: continue
      if not os.path.isdir(os.path.join(srcdir,'manualpages',llmansec)): os.mkdir(os.path.join(srcdir,'manualpages',llmansec))

      command = [doctext,
                 '-myst',
                 '-mpath',    os.path.join(srcdir,'manualpages',llmansec),
                 '-heading',  'SLEPc',
                 '-defn',     os.path.join(srcdir,'manualpages','doctext','myst.def'),
                 '-indexdir', '../'+llmansec,
                 '-index',    os.path.join(srcdir,'manualpages','manualpages.cit'),
                 '-locdir',   dir[len(slepc_dir)+1:]+'/',
                 '-Wargdesc', os.path.join(srcdir,'manualpages','doctext','doctextcommon.txt'),
                 file]
      sp = subprocess.run(command, cwd=dir, capture_output=True, encoding='UTF-8', check=True)
      if sp.stdout and sp.stdout.find('WARNING') > -1:
        print(sp.stdout)
        numberErrors = numberErrors + 1
      if sp.stderr and sp.stderr.find('WARNING') > -1:
        print(sp.stderr)
        numberErrors = numberErrors + 1
  return numberErrors


def processkhash(T, t, KeyType, ValType, text):
  '''Replaces T, t, KeyType, and ValType in text (from include/petsc/private/hashset.txt) with a set of supported values'''
  import re
  return re.sub('<ValType>',ValType,re.sub('<KeyType>',KeyType,re.sub('<t>',t,re.sub('<T>',T,text))))

def main(slepc_dir, srcdir, doctext):
  # generate the .md files for the manual pages from all the SLEPc source code
  try:
    os.unlink(os.path.join(srcdir,'manualpages','manualpages.cit'))
  except:
    pass
  numberErrors = 0
  for dirpath, dirnames, filenames in os.walk(os.path.join(slepc_dir),topdown=True):
    dirnames[:] = [d for d in dirnames if d not in ['tests', 'tutorials',
                                                    'doc', 'output',
                                                    'ftn-custom', 'ftn-auto',
                                                    'ftn-mod', 'binding',
                                                    'config', 'lib', '.git',
                                                    '.gitlab', 'share', 'systems'] and not d.startswith('arch')]
    numberErrors = numberErrors + processdir(slepc_dir,srcdir,dirpath,doctext)
  if numberErrors:
    raise RuntimeError('Stopping document build since errors were detected in generating manual pages')

  # generate list of all manual pages
  with open(os.path.join(srcdir,'manualpages','htmlmap'),mode='w') as map:
    with open(os.path.join(srcdir,'manualpages','manualpages.cit')) as cit:
      map.write(re.sub(r'man\+../','man+manualpages/',cit.read()))
    with open(os.path.join(srcdir,'manualpages','mpi.www.index')) as mpi:
      map.write(mpi.read())
