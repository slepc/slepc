#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#

import os, sys

class ArgDB:

  def __init__(self,argv):
    # standardize options
    for l in range(1,len(argv)):
      if argv[l].endswith('=yes'): argv[l] = argv[l].replace('=yes','=1')
      if argv[l].endswith('=no'): argv[l] = argv[l].replace('=no','=0')
    for l in range(1,len(argv)):
      name = argv[l]
      if name.startswith('--enable'):
        argv[l] = name.replace('--enable','--with')
        if name.find('=') == -1: argv[l] += '=1'
      elif name.startswith('--disable'):
        argv[l] = name.replace('--disable','--with')
        if name.find('=') == -1: argv[l] += '=0'
        elif name.endswith('=1'): argv[l] = argv[l].replace('=1','=0')
      elif name.startswith('--without'):
        argv[l] = name.replace('--without','--with')
        if name.find('=') == -1: argv[l] += '=0'
        elif name.endswith('=1'): argv[l] = argv[l].replace('=1','=0')
      elif name.startswith('--with'):
        if name.find('=') == -1: argv[l] += '=1'
      elif name.startswith('--download'):
        if name.endswith('=1'): argv[l] = name[:-2]
    self.argdb = argv[1:]
    self.useda = []

  def UsedArgs(self):
    return ' '.join(self.useda)

  def UsedArgsList(self):
    return self.useda

  def PopString(self,keyword):
    string = ''
    numhits = 0
    while True:
      found = 0
      for i, s in enumerate(self.argdb):
        if s.startswith('--'+keyword+'='):
          string = s.split('=',1)[1]
          found = 1
          numhits = numhits + 1
          self.useda.append(self.argdb[i])
          del self.argdb[i]
          break
      if not found:
        break
    return string,numhits

  def PopPath(self,keyword,exist=False):
    string = ''
    numhits = 0
    while True:
      found = 0
      for i, s in enumerate(self.argdb):
        if s.startswith('--'+keyword+'='):
          string = os.path.expanduser(s.split('=')[1].rstrip('/'))
          found = 1
          numhits = numhits + 1
          self.useda.append(self.argdb[i])
          del self.argdb[i]
          break
      if not found:
        break
    if exist and numhits and not os.path.exists(string):
      sys.exit('ERROR: The path given to --'+keyword+' does not exist: '+string)
    return string,numhits

  def PopUrl(self,keyword):
    value = False
    string = ''
    numhits = 0
    while True:
      found = 0
      for i, s in enumerate(self.argdb):
        if s.startswith('--'+keyword):
          value = not s.endswith('=0')
          try: string = s.split('=')[1]
          except IndexError: pass
          found = 1
          numhits = numhits + 1
          self.useda.append(self.argdb[i])
          del self.argdb[i]
          break
      if not found:
        break
    return string,value,numhits

  def PopBool(self,keyword):
    value = False
    numhits = 0
    while True:
      found = 0
      for i, s in enumerate(self.argdb):
        if s == '--'+keyword or s.startswith('--'+keyword+'='):
          value = not s.endswith('=0')
          found = 1
          numhits = numhits + 1
          self.useda.append(self.argdb[i])
          del self.argdb[i]
          break
      if not found:
        break
    return value,numhits

  def PopHelp(self):
    value = False
    numhits = 0
    while True:
      found = 0
      for i, s in enumerate(self.argdb):
        s = s.strip()
        if s in ['--help', '-h', '-?']:
          value = True
          found = 1
          numhits = numhits + 1
          del self.argdb[i]
          break
      if not found:
        break
    return value

  def ErrorPetscOptions(self):
    petscopts = []
    strings = ['with-precision','with-scalar-type','with-clanguage','with-cc','with-cxx','with-fc','with-blaslapack-lib','with-scalapack-dir','with-elemental-dir']
    for s in strings:
      value,found = self.PopString(s)
      if found: petscopts.append(s)
    bools = ['with-fortran-bindings','with-64-bit-indices','with-shared-libraries','with-single-library','with-debugging','with-cuda','with-hip','with-mpi','with-openmp']
    for s in bools:
      value,found = self.PopBool(s)
      if found: petscopts.append(s)
    urls = ['download-f2cblaslapack','download-mumps','download-scalapack','download-elemental']
    for s in urls:
      url,flag,found = self.PopUrl(s)
      if found: petscopts.append(s)
    if petscopts:
      sys.exit('ERROR: The following options belong to PETSc configure: '+', '.join(petscopts)+'\nUse -h for help')

  def ErrorIfNotEmpty(self):
    if self.argdb:
      sys.exit('ERROR: Invalid arguments '+' '.join(self.argdb)+'\nUse -h for help')

