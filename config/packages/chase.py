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

class Chase(package.Package):

  def __init__(self,argdb,log):
    package.Package.__init__(self,argdb,log)
    self.packagename    = 'chase'
    self.packagetype    = 'cmake'
    self.installable    = True
    self.downloadable   = True
    self.version        = '1.5.0'
    obj = self.version if hasattr(self,'version') else self.gitcommit
    self.url            = 'https://github.com/ChASE-library/ChASE/archive/'+('v'+obj if hasattr(self,'version') else obj)+'.tar.gz'
    self.archive        = 'ChASE-'+obj+'.tar.gz'
    self.supportssingle = True
    self.hasheaders     = True
    self.ProcessArgs(argdb)

  def Precondition(self,slepc,petsc):
    pkg = self.packagename.upper()
    if not 'scalapack' in petsc.packages:
      self.log.Exit(pkg+' requires PETSc to be configured with ScaLAPACK')
    if hasattr(self,'download') and self.download:
      if not hasattr(petsc,'cmake'):
        self.log.Exit(pkg+' requires CMake for building')
      if petsc.maxcxxdialect == '':
        self.log.Exit('Downloading '+pkg+' requires a functioning C++ compiler')
    package.Package.Precondition(self,slepc,petsc)

  def SampleCode(self,petsc):
    code =  '#include <mpi.h>\n'
    code += 'void zchase_(int* deg, double* tol, char* mode, char* opt, char *qr);\n'
    code += 'int main(int argc, char** argv) {\n'
    code += '  int deg;\n'
    code += '  char mode,opt,qr;\n'
    code += '  double tol;\n'
    code += '  MPI_Init(&argc, &argv);\n'
    code += '  zchase_(&deg, &tol, &mode, &opt, &qr);\n'
    code += '  MPI_Finalize();\n'
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
      dirs = self.GenerateGuesses('chase',archdir) + self.GenerateGuesses('chase',archdir,'lib64')
      incdirs = self.GenerateGuesses('chase',archdir,'include')

    libs = [self.packagelibs] if self.packagelibs else ['-lchase_c']
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
          slepcconf.write('#define SLEPC_HAVE_CHASE 1\n')
          slepcvars.write('CHASE_LIB = ' + self.libflags + '\n')
          slepcvars.write('CHASE_INCLUDE = ' + self.includeflags + '\n')
          self.havepackage = True
          self.packageflags = ' '.join(l+f)
          return

    self.log.Exit('Unable to link with ChASE library in directories'+' '.join(dirs)+' with libraries and link flags '+' '.join(libs))


  def DownloadAndInstall(self,slepcconf,slepcvars,slepc,petsc,archdir,prefixdir):
    externdir = slepc.GetExternalPackagesDir(archdir)
    builddir  = self.Download(externdir,slepc.downloaddir)
    incdir,libdir = slepc.CreatePrefixDirs(prefixdir)

    builddir = slepc.CreateDir(builddir,'build')
    confopt = ['-DCMAKE_INSTALL_PREFIX='+prefixdir, '-DCMAKE_INSTALL_NAME_DIR:STRING="'+os.path.join(prefixdir,'lib')+'"', '-DCMAKE_INSTALL_LIBDIR:STRING="lib"', '-DCMAKE_C_COMPILER="'+petsc.cc+'"', '-DCMAKE_C_FLAGS:STRING="'+petsc.getCFlags()+'"', '-DCMAKE_CXX_COMPILER="'+petsc.cxx+'"', '-DCMAKE_CXX_FLAGS:STRING="'+petsc.getCXXFlags()+'"', '-DCMAKE_Fortran_COMPILER="'+petsc.fc+'"', '-DCMAKE_Fortran_FLAGS:STRING="'+petsc.getFFlags()+'"', '-DSCALAPACK_LIBRARIES="'+petsc.scalapack_lib+'"', '-DBLAS_LIBRARIES="'+petsc.blaslapack_lib+'"']
    confopt.append('-DCMAKE_BUILD_TYPE='+('Debug' if petsc.debug else 'Release'))
    if petsc.buildsharedlib:
      confopt = confopt + ['-DBUILD_SHARED_LIBS=ON', '-DCMAKE_INSTALL_RPATH:PATH='+os.path.join(prefixdir,'lib')]
    else:
      confopt.append('-DBUILD_SHARED_LIBS=OFF')
    if 'MSYSTEM' in os.environ:
      confopt.append('-G "MSYS Makefiles"')
    (result,output) = self.RunCommand('cd '+builddir+' && '+petsc.cmake+' '+' '.join(confopt)+' '+self.buildflags+' .. && '+petsc.make+' -j'+petsc.make_np+' && '+petsc.make+' install')

    if result:
      self.log.Exit('Installation of ChASE failed')

    # Check build
    code = self.SampleCode(petsc)
    altlibdir = os.path.join(prefixdir,'lib64')
    for ldir in [libdir,altlibdir]:
      if petsc.buildsharedlib:
        l = self.slflag + ldir + ' -L' + ldir + ' -lchase_c'
      else:
        l = '-L' + ldir + ' -lchase_c'
      f = '-I' + incdir
      (result,output) = self.Link([],[],l+' '+f,code,f,petsc.language)
      if result: break

    if not result:
      self.log.Exit('Unable to link with downloaded ChASE')

    # Write configuration files
    self.libflags = l
    self.includeflags = f
    slepcconf.write('#define SLEPC_HAVE_CHASE 1\n')
    slepcvars.write('CHASE_LIB = ' + self.libflags + '\n')
    slepcvars.write('CHASE_INCLUDE = ' + self.includeflags + '\n')

    self.havepackage = True
    self.packageflags = l+' '+f

