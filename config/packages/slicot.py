#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#

import os, log, package

class Slicot(package.Package):

  def __init__(self,argdb,log):
    package.Package.__init__(self,argdb,log)
    self.packagename    = 'slicot'
    self.packagetype    = 'cmake'
    self.installable    = True
    self.downloadable   = True
    #self.gitcommit      = 'a037f7eb76134d45e7d222b7f017d5cbd16eb731'
    self.version        = '5.9.1'
    obj = self.version if hasattr(self,'version') else self.gitcommit
    self.url            = 'https://github.com/SLICOT/SLICOT-Reference/archive/'+('v'+obj if hasattr(self,'version') else obj)+'.tar.gz'
    self.archive        = 'slicot-'+obj+'.tar.gz'
    self.supportsscalar = ['real']
    self.fortran        = True
    self.ProcessArgs(argdb)

  def Check(self,slepcconf,slepcvars,petsc,archdir):
    functions = ['sb03od','sb03md']
    libs = self.packagelibs if self.packagelibs else [['-lslicot']]

    if self.packagedir:
      if os.path.isdir(os.path.join(os.sep,'usr','lib64')):
        dirs = ['',os.path.join(self.packagedir,'lib64'),self.packagedir,os.path.join(self.packagedir,'lib')]
      else:
        dirs = ['',os.path.join(self.packagedir,'lib'),self.packagedir,os.path.join(self.packagedir,'lib64')]
    else:
      dirs = self.GenerateGuesses('slicot',archdir)

    self.FortranLib(slepcconf,slepcvars,dirs,libs,functions)

  def DownloadAndInstall(self,slepcconf,slepcvars,slepc,petsc,archdir,prefixdir):
    externdir = slepc.GetExternalPackagesDir(archdir)
    builddir  = self.Download(externdir,slepc.downloaddir)

    # Build with cmake
    builddir = slepc.CreateDir(builddir,'build')
    confopt = ['-DCMAKE_INSTALL_PREFIX='+prefixdir, '-DCMAKE_INSTALL_NAME_DIR:STRING="'+os.path.join(prefixdir,'lib')+'"', '-DCMAKE_INSTALL_LIBDIR:STRING="lib"', '-DCMAKE_C_COMPILER="'+petsc.cc+'"', '-DCMAKE_C_FLAGS:STRING="'+petsc.getCFlags()+'"', '-DCMAKE_Fortran_COMPILER="'+petsc.fc+'"', '-DCMAKE_Fortran_FLAGS:STRING="'+petsc.getFFlags()+'"', '-DSLICOT_TESTING=OFF']
    confopt.append('-DCMAKE_BUILD_TYPE='+('Debug' if petsc.debug else 'Release'))
    if petsc.buildsharedlib:
      confopt = confopt + ['-DBUILD_SHARED_LIBS=ON', '-DCMAKE_INSTALL_RPATH:PATH='+os.path.join(prefixdir,'lib')]
    else:
      confopt.append('-DBUILD_SHARED_LIBS=OFF')
    if petsc.ind64:
      confopt.append('-DSLICOT_INTEGER8=ON')
    if 'MSYSTEM' in os.environ:
      confopt.append('-G "MSYS Makefiles"')
    (result,output) = self.RunCommand('cd '+builddir+' && '+petsc.cmake+' '+' '.join(confopt)+' '+self.buildflags+' .. && '+petsc.make+' -j'+petsc.make_np+' && '+petsc.make+' install')

    # Check build
    functions = ['sb03od']
    libs = [['-lslicot']]
    dirs = [os.path.join(prefixdir,'lib'),os.path.join(prefixdir,'lib64')]
    self.FortranLib(slepcconf,slepcvars,dirs,libs,functions)
