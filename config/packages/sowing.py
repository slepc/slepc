#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#

from __future__ import print_function
import os, sys, log, package

class Sowing(package.Package):

  def __init__(self,argdb,log):
    package.Package.__init__(self,argdb,log)
    self.packagename  = 'sowing'
    self.downloadable = True
    #self.gitcommit    = '3a410fa51a7bb531676f16deb8bc0c1ded8293c3'
    self.version      = '1.1.26.12'
    obj = self.version if hasattr(self,'version') else self.gitcommit
    self.url          = 'https://bitbucket.org/petsc/pkg-sowing/get/'+('v'+obj if hasattr(self,'version') else obj)+'.tar.gz'
    self.archive      = 'sowing-'+obj+'.tar.gz'
    self.ProcessArgs(argdb)

  def ShowHelp(self):
    wd = package.Package.wd
    print('  --download-sowing[=<fname>]'.ljust(wd)+': Download and install SOWING (developers, to build documentation)')

  def DownloadAndInstall(self,slepcconf,slepcvars,slepc,petsc,archdir,prefixdir):
    externdir = slepc.GetExternalPackagesDir(archdir)
    builddir  = self.Download(externdir,slepc.downloaddir)

    # Configure, build and install package
    (result,output) = self.RunCommand('cd '+builddir+'&& ./configure --prefix='+prefixdir+'&&'+petsc.make+'&&'+petsc.make+' install')

    slepcvars.write('DOCTEXT = ' + os.path.join(prefixdir,'bin','doctext') + '\n')
    slepcvars.write('MAPNAMES = ' + os.path.join(prefixdir,'bin','mapnames') + '\n')
    self.havepackage = True

