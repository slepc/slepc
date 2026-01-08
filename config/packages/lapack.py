#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#

import log, package

class Lapack(package.Package):

  def __init__(self,argdb,log):
    package.Package.__init__(self,argdb,log)
    self.packagename = 'lapack'

  def ShowHelp(self):
    print('BLAS/LAPACK:')
    print('  [Not optional, always checked]')

  def ShowInfo(self):
    if hasattr(self,'missing'):
      self.log.Warn('Missing LAPACK functions: %s.\nSome SLEPc functionality will not be available.\nPlease reconfigure and recompile PETSc with a full LAPACK implementation'%(' '.join(self.missing)))

  def Process(self,slepcconf,slepcvars,slepcrules,slepc,petsc,archdir=''):
    self.make = petsc.make
    self.mangling = petsc.blaslapackmangling
    if petsc.buildsharedlib:
      self.slflag = petsc.cc_linker_slflag
    self.log.NewSection('Checking LAPACK library...')
    self.Check(slepcconf,slepcvars,petsc)

  def LinkBlasLapack(self,functions,callbacks,flags,petsc):
    code = ''
    for f in functions:
      if petsc.language == 'c++':
        code += 'extern "C" void ' + f + '();\n'
      else:
        code += 'extern void ' + f + '();\n'
    code += 'int main() {\n'
    for f in functions:
      code += f + '();\n'
    code += 'return 0;\n}\n'
    (result, output) = self.Link(functions,callbacks,flags,code)
    return result

  def Mangle(self,i):
    if self.mangling == 'underscore':
      f = i + '_'
    elif self.mangling == 'caps':
      f = i.upper()
    else:
      f = i
    return f

  def Check(self,slepcconf,slepcvars,petsc):

    # LAPACK standard functions
    l = ['laev2','gehrd','lanhs','trexc','trevc','tgexc','tgevc','stedc','hsein','larfg','larf','lascl','trsyl','bdsqr','gebrd']

    # LAPACK functions with different real and complex names
    if petsc.scalar == 'real':
      l += ['orghr','sytrd','orgtr','ormbr']
      if petsc.precision == 'single':
        prefix = 's'
      elif petsc.precision == '__float128':
        prefix = 'q'
      else:
        prefix = 'd'
    else:
      l += ['unghr','hetrd','ungtr','unmbr']
      if petsc.precision == 'single':
        prefix = 'c'
      elif petsc.precision == '__float128':
        prefix = 'w'
      else:
        prefix = 'z'

    # add prefix to LAPACK names
    functions = [prefix + i for i in l]

    # in this case, the real name represents both versions
    namesubst = {'unghr':'orghr','hetrd':'sytrd','ungtr':'orgtr','unmbr':'ormbr'}

    # LAPACK functions which are always used in real version
    l = ['stevr','bdsdc','lag2','lasv2','lartg','laed4','lamrg']
    if petsc.precision == 'single':
      rprefix = 's'
    elif petsc.precision == '__float128':
      rprefix = 'q'
    else:
      rprefix = 'd'
    for i in l:
      functions.append(rprefix + i)

    # check for all functions at once
    allf = [self.Mangle(i) for i in functions]

    self.log.write('=== Checking all LAPACK functions...')
    if not self.LinkBlasLapack(allf,[],'',petsc):
      # check functions one by one
      self.missing = []
      for i in functions:
        self.log.write('=== Checking LAPACK '+i+' function...')
        if not self.LinkBlasLapack([self.Mangle(i)],[],'',petsc):
          self.missing.append(i)
          # some complex functions are represented by their real names
          if i[1:] in namesubst:
            nf = namesubst[i[1:]]
          else:
            nf = i[1:]
          slepcconf.write('#define SLEPC_MISSING_LAPACK_' + nf.upper() + ' 1\n')

    # check for pairs of advanced and basic subroutines
    # if advanced is missing (only available in recent Lapack) then check the basic one
    for pair in [['ggsvd3','ggsvd'],['ggev3','ggev'],['gges3','']]:   # gges is checked in PETSc
      for f in pair:
        if not f: continue
        i = prefix + f
        self.log.write('=== Checking LAPACK '+i+' function...')
        if not self.LinkBlasLapack([self.Mangle(i)],[],'',petsc):
          if f == pair[1]: # do not warn if advanced routine is missing
            if hasattr(self,'missing'): self.missing.append(i)
            else: self.missing = [i]
          nf = i[1:]
          slepcconf.write('#define SLEPC_MISSING_LAPACK_' + nf.upper() + ' 1\n')
        else: break
