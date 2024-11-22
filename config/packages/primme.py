#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#

import os, shutil, log, package

class Primme(package.Package):

  def __init__(self,argdb,log):
    package.Package.__init__(self,argdb,log)
    self.packagename    = 'primme'
    self.packagetype    = 'source_c'
    self.installable    = True
    self.downloadable   = True
    self.version        = '3.2.3'
    self.url            = 'https://github.com/primme/primme/archive/v'+self.version+'.tar.gz'
    self.archive        = 'primme-'+self.version+'.tar.gz'
    self.supportssingle = True
    self.supports64bint = True
    self.hasheaders     = True
    self.ProcessArgs(argdb)

  def SampleCode(self,petsc):
    if petsc.scalar == 'real':
      if petsc.precision == 'single':
        function = 'sprimme'
        rdtype = 'float'
      else:
        function = 'dprimme'
        rdtype = 'double'
      cdtype = rdtype
    else:
      if petsc.precision == 'single':
        function = 'cprimme'
        rdtype = 'float'
        cdtype = 'PRIMME_COMPLEX_FLOAT'
      else:
        function = 'zprimme'
        rdtype = 'double'
        cdtype = 'PRIMME_COMPLEX_DOUBLE'

    code = '#include "primme.h"\n'
    code += 'int main() {\n'
    code += '  ' + rdtype + ' *a=NULL,*c=NULL;\n'
    code += '  ' + cdtype + ' *b=NULL;\n'
    code += '  primme_params primme;\n'
    code += '  primme_initialize(&primme);\n'
    code += '  primme_set_method(PRIMME_DYNAMIC,&primme);\n'
    code += '  ' + function + '(a,b,c,&primme);\n'
    code += '  primme_free(&primme);\n'
    code += '  return 0;\n}\n'
    return code


  def Check(self,slepcconf,slepcvars,petsc,archdir):
    code = self.SampleCode(petsc)
    if self.packagedir:
      if os.path.isdir(os.path.join(os.sep,'usr','lib64')):
        dirs = ['',os.path.join(self.packagedir,'lib64'),self.packagedir,os.path.join(self.packagedir,'lib')]
      else:
        dirs = ['',os.path.join(self.packagedir,'lib'),self.packagedir,os.path.join(self.packagedir,'lib64')]
      incdirs = ['',os.path.join(self.packagedir,'include'),self.packagedir]
    else:
      dirs = self.GenerateGuesses('Primme',archdir) + self.GenerateGuesses('Primme',archdir,'lib64')
      incdirs = self.GenerateGuesses('Primme',archdir,'include')

    libs = [self.packagelibs] if self.packagelibs else ['-lprimme']
    includes = [self.packageincludes] if self.packageincludes else []

    for d in dirs:
      for i in incdirs:
        if d:
          if petsc.buildsharedlib:
            l = [self.slflag + d] + ['-L' + d] + libs
          else:
            l = ['-L' + d] + libs
          f = (['-I' + i] if i else [])
        else:
          l = libs
          f = []
        (result, output) = self.Link([],[],' '.join(l+f),code,' '.join(f),petsc.language)
        if result:
          self.libflags = ' '.join(l)
          self.includeflags = ' '.join(f)
          slepcconf.write('#define SLEPC_HAVE_PRIMME 1\n')
          slepcvars.write('PRIMME_LIB = ' + self.libflags + '\n')
          slepcvars.write('PRIMME_INCLUDE = ' + self.includeflags + '\n')
          self.havepackage = True
          self.packageflags = ' '.join(l+f)
          self.location = includes[0] if self.packageincludes else i
          if not self.location:
            self.location = self.DefaultIncludePath(petsc,'primme.h')
          return

    self.log.Exit('Unable to link with PRIMME library in directories'+' '.join(dirs)+' with libraries and link flags '+' '.join(libs)+' [NOTE: make sure PRIMME version is 2.0 at least]')


  def DownloadAndInstall(self,slepcconf,slepcvars,slepc,petsc,archdir,prefixdir):
    externdir = slepc.GetExternalPackagesDir(archdir)
    builddir  = self.Download(externdir,slepc.downloaddir)

    # Makefile
    cont =  'export LIBRARY       = libprimme.'+petsc.ar_lib_suffix+'\n'
    cont += 'export SOLIBRARY     = libprimme.'+petsc.sl_linker_suffix+'\n'
    cont += 'export SONAMELIBRARY = libprimme.'+petsc.sl_linker_suffix+'.'+self.version+'\n'
    cont += 'export CC            = '+petsc.cc+'\n'
    if hasattr(petsc,'fc'):
      cont += 'export F77           = '+petsc.fc+'\n'
    cont += 'export DEFINES       = '
    if petsc.blaslapackmangling == 'underscore':
      cont += '-DF77UNDERSCORE '
    if petsc.blaslapackint64:
      cont += '-DPRIMME_BLASINT_SIZE=64'
    cont += '\n'
    cont += 'export INCLUDE       = \n'
    cont += 'export CFLAGS        = '+petsc.getCFlags()+' '+self.buildflags+'\n'
    cont += 'export RANLIB        = '+petsc.ranlib+'\n'
    cont += 'export PREFIX        = '+prefixdir+'\n'
    cont += 'include makefile\n'
    self.WriteMakefile('mymake_flags',builddir,cont)

    # Build package
    target = ' install' if petsc.buildsharedlib else ' lib'
    mymake = petsc.make + ' -j' + petsc.make_np + ' -f mymake_flags '
    (result, output) = self.RunCommand('cd '+builddir+'&&'+mymake+' clean && '+mymake+target)
    if result:
      self.log.Exit('Installation of PRIMME failed')

    # Move files
    incdir,libdir = slepc.CreatePrefixDirs(prefixdir)
    if not petsc.buildsharedlib:
      os.rename(os.path.join(builddir,'lib','libprimme.'+petsc.ar_lib_suffix),os.path.join(libdir,'libprimme.'+petsc.ar_lib_suffix))
      for root, dirs, files in os.walk(os.path.join(builddir,'include')):
        for name in files:
          shutil.copyfile(os.path.join(builddir,'include',name),os.path.join(incdir,name))

    if petsc.buildsharedlib:
      l = self.slflag + libdir + ' -L' + libdir + ' -lprimme'
    else:
      l = '-L' + libdir + ' -lprimme'
    f = '-I' + incdir

    # Check build
    code = self.SampleCode(petsc)
    (result, output) = self.Link([],[],l+' '+f,code,f,petsc.language)
    if not result:
      self.log.Exit('Unable to link with downloaded PRIMME')

    # Write configuration files
    self.libflags = l
    self.includeflags = f
    slepcconf.write('#define SLEPC_HAVE_PRIMME 1\n')
    slepcvars.write('PRIMME_LIB = ' + self.libflags + '\n')
    slepcvars.write('PRIMME_INCLUDE = ' + self.includeflags + '\n')

    self.location = incdir
    self.havepackage = True
    self.packageflags = l+' '+f


  def LoadVersion(self,slepcconf):
    try:
      with open(os.path.join(self.location,'primme.h')) as f:
        for l in f.readlines():
          l = l.split()
          if len(l) == 3:
            if l[1] == 'PRIMME_VERSION_MAJOR':
              major = l[2]
            elif l[1] == 'PRIMME_VERSION_MINOR':
              minor = l[2]
      self.iversion = major + '.' + minor
      if major=='3':
        slepcconf.write('#define SLEPC_HAVE_PRIMME3 1\n')
    except Exception as e:
      self.log.write('Error while determining version of PRIMME:\n'+str(e))

