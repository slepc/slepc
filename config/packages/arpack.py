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

class Arpack(package.Package):

  def __init__(self,argdb,log):
    package.Package.__init__(self,argdb,log)
    self.packagename    = 'arpack'
    self.packagetype    = 'cmake'
    self.installable    = True
    self.downloadable   = True
    self.version        = '3.9.1'
    #self.gitcommit      = '5131f792f289c4e63b4cb1f56003e59507910132'
    obj = self.version if hasattr(self,'version') else self.gitcommit
    self.url            = 'https://github.com/opencollab/arpack-ng/archive/'+obj+'.tar.gz'
    self.archive        = 'arpack-ng-'+obj+'.tar.gz'
    self.supports64bint = True
    self.fortran        = True
    self.hasheaders     = True   # the option --with-arpack-include=... is simply ignored
    self.supportsprecis.append('single')
    self.ProcessArgs(argdb)

  def Functions(self,petsc):
    if petsc.mpiuni or petsc.msmpi:
      if petsc.scalar == 'real':
        if petsc.precision == 'single':
          functions = ['snaupd','sneupd','ssaupd','sseupd']
        else:
          functions = ['dnaupd','dneupd','dsaupd','dseupd']
      else:
        if petsc.precision == 'single':
          functions = ['cnaupd','cneupd']
        else:
          functions = ['znaupd','zneupd']
    else:
      if petsc.scalar == 'real':
        if petsc.precision == 'single':
          functions = ['psnaupd','psneupd','pssaupd','psseupd']
        else:
          functions = ['pdnaupd','pdneupd','pdsaupd','pdseupd']
      else:
        if petsc.precision == 'single':
          functions = ['pcnaupd','pcneupd']
        else:
          functions = ['pznaupd','pzneupd']
    return functions

  def Check(self,slepcconf,slepcvars,petsc,archdir):
    functions = self.Functions(petsc)
    if self.packagelibs:
      libs = self.packagelibs
    else:
      if petsc.mpiuni or petsc.msmpi:
        libs = [['-larpack'],['-larpack_LINUX'],['-larpack_SUN4']]
      else:
        libs = [['-lparpack','-larpack'],['-lparpack_MPI','-larpack'],['-lparpack_MPI-LINUX','-larpack_LINUX'],['-lparpack_MPI-SUN4','-larpack_SUN4']]

    if self.packagedir:
      if os.path.isdir(os.path.join(os.sep,'usr','lib64')):
        dirs = ['',os.path.join(self.packagedir,'lib64'),self.packagedir,os.path.join(self.packagedir,'lib')]
      else:
        dirs = ['',os.path.join(self.packagedir,'lib'),self.packagedir,os.path.join(self.packagedir,'lib64')]
    else:
      dirs = self.GenerateGuesses('Arpack',archdir) + self.GenerateGuesses('Arpack',archdir,'lib64')
    self.FortranLib(slepcconf,slepcvars,dirs,libs,functions)

  def DownloadAndInstall(self,slepcconf,slepcvars,slepc,petsc,archdir,prefixdir):
    externdir = slepc.GetExternalPackagesDir(archdir)
    builddir  = self.Download(externdir,slepc.downloaddir)

    # Check user options
    if petsc.ind64:
      if not petsc.mpiuni:
        self.log.Exit('Parallel ARPACK does not support 64-bit integers')
      if not petsc.blaslapackint64:
        self.log.Exit('To install ARPACK with 64-bit integers you also need a BLAS with 64-bit integers')

    if hasattr(petsc,'cmake'): # Build with cmake
      builddir = slepc.CreateDir(builddir,'build')
      confopt = ['-DCMAKE_INSTALL_PREFIX='+prefixdir, '-DCMAKE_INSTALL_NAME_DIR:STRING="'+os.path.join(prefixdir,'lib')+'"', '-DCMAKE_INSTALL_LIBDIR:STRING="lib"', '-DCMAKE_C_COMPILER="'+petsc.cc+'"', '-DCMAKE_C_FLAGS:STRING="'+petsc.getCFlags()+'"', '-DCMAKE_Fortran_COMPILER="'+petsc.fc+'"', '-DCMAKE_Fortran_FLAGS:STRING="'+petsc.getFFlags()+'"', '-DBLAS_LIBRARIES="'+petsc.blaslapack_lib+'"']
      if not petsc.mpiuni and not petsc.msmpi:
        confopt = confopt + ['-DMPI=ON', '-DMPI_C_COMPILER="'+petsc.cc+'"', '-DMPI_Fortran_COMPILER="'+petsc.fc+'"']
      confopt.append('-DCMAKE_BUILD_TYPE='+('Debug' if petsc.debug else 'Release'))
      if petsc.buildsharedlib:
        confopt = confopt + ['-DBUILD_SHARED_LIBS=ON', '-DCMAKE_INSTALL_RPATH:PATH='+os.path.join(prefixdir,'lib')]
      else:
        confopt.append('-DBUILD_SHARED_LIBS=OFF')
      if petsc.ind64:
        confopt.append('-DINTERFACE64=1')
      if 'MSYSTEM' in os.environ:
        confopt.append('-G "MSYS Makefiles"')
      (result,output) = self.RunCommand('cd '+builddir+' && '+petsc.cmake+' '+' '.join(confopt)+' '+self.buildflags+' .. && '+petsc.make+' -j'+petsc.make_np+' && '+petsc.make+' install')

    else: # Build with autoreconf
      (result,output) = self.RunCommand('autoreconf --help')
      if result:
        self.log.Exit('--download-arpack requires that the command autoreconf is available on your PATH, or alternatively that PETSc has been configured with CMake')
      if self.buildflags:
        self.log.Exit('You specified --download-arpack-cmake-arguments but ARPACK will be built with autoreconf because PETSc was not configured with CMake')
      confopt = ['--prefix='+prefixdir, '--libdir='+os.path.join(prefixdir,'lib'), 'CC="'+petsc.cc+'"', 'CFLAGS="'+petsc.getCFlags()+'"', 'F77="'+petsc.fc+'"', 'FFLAGS="'+petsc.getFFlags()+'"', 'FC="'+petsc.fc+'"', 'FCFLAGS="'+petsc.getFFlags()+'"', 'LIBS="'+petsc.blaslapack_lib+'"']
      if not petsc.mpiuni and not petsc.msmpi:
        confopt = confopt + ['--enable-mpi MPICC="'+petsc.cc+'"', 'MPIF77="'+petsc.fc+'"', 'MPIFC="'+petsc.fc+'"']
      if not petsc.buildsharedlib:
        confopt.append('--disable-shared')
      if petsc.ind64:
        confopt.append('INTERFACE64=1')
      (result,output) = self.RunCommand('cd '+builddir+'&& sh bootstrap && ./configure '+' '.join(confopt)+' && '+petsc.make+' -j'+petsc.make_np+' && '+petsc.make+' install')

    if result:
      self.log.Exit('Installation of ARPACK failed')

    # Check build
    functions = self.Functions(petsc)
    if petsc.mpiuni or petsc.msmpi:
      libs = [['-larpack']]
    else:
      libs = [['-lparpack','-larpack']]
    dirs = [os.path.join(prefixdir,'lib'),os.path.join(prefixdir,'lib64')]
    self.FortranLib(slepcconf,slepcvars,dirs,libs,functions)
