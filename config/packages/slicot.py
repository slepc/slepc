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
    self.installable    = True
    self.downloadable   = True
    #self.gitcommit      = 'a037f7eb76134d45e7d222b7f017d5cbd16eb731'
    self.version        = '5.9'
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
    libname = 'libslicot.a'

    # Makefile
    cont  = 'FORTRAN   = '+petsc.fc+'\n'
    cont += 'OPTS      = '+petsc.getFFlags()+'\n'
    cont += 'ARCH      = '+petsc.ar+'\n'
    cont += 'ARCHFLAGS = '+petsc.ar_flags+'\n'
    cont += 'SLICOTLIB = ../'+libname+'\n'
    cont += 'LPKAUXLIB = ../'+libname+'\n'  # TODO: use a separate library for this
    self.WriteMakefile('make_Unix.inc',builddir,cont)

    # Build package
    target = 'lib'
    (result,output) = self.RunCommand('cd '+builddir+' && '+petsc.make+' -f makefile_Unix cleanlib && '+petsc.make+' -f makefile_Unix -j'+petsc.make_np+' '+target)
    if result:
      self.log.Exit('Installation of SLICOT failed')

    # Move files
    incdir,libdir = slepc.CreatePrefixDirs(prefixdir)
    os.rename(os.path.join(builddir,libname),os.path.join(libdir,libname))

    # Check build
    functions = ['sb03od']
    libs = [['-lslicot']]
    dirs = [libdir]
    self.FortranLib(slepcconf,slepcvars,dirs,libs,functions)
