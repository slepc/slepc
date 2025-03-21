#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#

import os,shutil,log,package
from urllib.request import urlretrieve

class HPDDM(package.Package):

  def __init__(self,argdb,log):
    package.Package.__init__(self,argdb,log)
    self.packagename    = 'hpddm'
    self.downloadable   = True
    self.gitcommit      = 'f940190ed791984e5ca501a58eb8e6d112ac231d'
    # self.version        = '2.3.1'
    obj = self.version if hasattr(self,'version') else self.gitcommit
    self.url            = 'https://github.com/hpddm/hpddm/archive/'+('v'+obj if hasattr(self,'version') else obj)+'.tar.gz'
    self.archive        = 'hpddm-'+obj+'.tar.gz'
    self.supportssingle = True
    self.supports64bint = True
    self.ProcessArgs(argdb)

  def Precondition(self,slepc,petsc):
    pkg = self.packagename.upper()
    if petsc.maxcxxdialect == '':
      self.log.Exit(pkg+' requires a functioning C++ compiler')
    if not petsc.buildsharedlib:
      self.log.Exit(pkg+' requires a shared library build')
    if 'slepc' in petsc.packages:
      self.log.Exit(pkg+' requires PETSc to be built without SLEPc')
    package.Package.Precondition(self,slepc,petsc)

  def DownloadAndInstall(self,slepcconf,slepcvars,slepc,petsc,archdir,prefixdir):
    externdir = slepc.GetExternalPackagesDir(archdir)
    builddir  = self.Download(externdir,slepc.downloaddir)
    incdir,libdir = slepc.CreatePrefixDirs(prefixdir)
    cont  = 'include '+os.path.join(petsc.dir,petsc.arch,'lib','petsc','conf','petscvariables')+'\n'
    cont += 'soname:\n'
    cont += '\t@echo $(call SONAME_FUNCTION,'+os.path.join(libdir,'libhpddm_petsc')+',0)\n'
    cont += 'sl_linker:\n'
    cont += '\t@echo $(call SL_LINKER_FUNCTION,'+os.path.join(libdir,'libhpddm_petsc')+',0,0)\n'
    self.WriteMakefile('SONAME_SL_LINKER',builddir,cont)
    d = os.path.join(petsc.dir,petsc.arch,'lib')
    l = self.slflag+d+' -L'+d+' -lpetsc'
    d = libdir
    cmd = petsc.cxx+' '+petsc.getCXXFlags()+' -I'+os.path.join('.','include')+' -I'+os.path.join(petsc.dir,petsc.arch,'include')+' -I'+os.path.join(slepc.dir,'include')+' -I'+os.path.join(archdir,'include')+' -I'+os.path.join(petsc.dir,'include')+' -DPETSC_HAVE_SLEPC=1 -DSLEPC_LIB_DIR="'+d+'"'
    line = 'cd '+builddir+' && '+cmd+' '
    if self.packagename not in petsc.packages: # KSPHPDDM and PCHPDDM sources have not been compiled by PETSc which was configured without --download-hpddm, so do it here
      line = line+os.path.join('interface','petsc','ksp','hpddm.cxx')+' -c -o '+os.path.join('interface','ksphpddm.o')+' && '+cmd+' '+os.path.join('interface','petsc','pc','pchpddm.cxx')+' -c -o '+os.path.join('interface','pchpddm.o')+' && '+cmd+' '
    line = line+os.path.join('interface','hpddm_petsc.cpp')+' -c -o '+os.path.join('interface','hpddm_petsc.o')
    (result,output) = self.RunCommand(line)
    if result:
      self.log.Exit('Compilation of HPDDM failed')
    (result,output) = self.RunCommand('cd '+builddir+' && make -f SONAME_SL_LINKER soname && make -f SONAME_SL_LINKER sl_linker')
    if result:
      self.log.Exit('Calling PETSc SONAME_FUNCTION or SL_LINKER_FUNCTION failed')
    lines = output.splitlines()
    line = 'cd '+builddir+' && '+petsc.cxx+' '+petsc.getCXXFlags()+' '+os.path.join('interface','hpddm_petsc.o')
    if self.packagename not in petsc.packages: # link KSPHPDDM and PCHPDDM objects in libhpddm_petsc since it was not done when building libpetsc
      line = line+' '+os.path.join('interface','pchpddm.o')+' '+os.path.join('interface','ksphpddm.o')
    line = line+' -o '+lines[0]+' '+lines[1]+' '+l+' && ln -sf '+lines[0]+' '+os.path.join(d,'libhpddm_petsc.'+petsc.sl_linker_suffix)
    (result,output) = self.RunCommand(line)
    if result:
      self.log.Exit('Installation of HPDDM failed')
    for root,dirs,files in os.walk(os.path.join(builddir,'include')):
      for name in files:
        shutil.copyfile(os.path.join(builddir,'include',name),os.path.join(incdir,name))
    l = self.slflag+d+' -L'+d+' -lhpddm_petsc'
    f = '-I'+incdir
    # Write configuration files
    self.libflags = l
    self.includeflags = f
    slepcconf.write('#define SLEPC_HAVE_HPDDM 1\n')
    slepcvars.write('HPDDM_LIB = '+self.libflags+'\n')
    slepcvars.write('HPDDM_INCLUDE = '+self.includeflags+'\n')
    self.packageflags = l+' '+f
    self.havepackage = True
