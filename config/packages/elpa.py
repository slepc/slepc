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

class Elpa(package.Package):

  def __init__(self,argdb,log):
    package.Package.__init__(self,argdb,log)
    self.packagename    = 'elpa'
    self.packagetype    = 'gnu'
    self.installable    = True
    self.downloadable   = True
    self.version        = '2025.01.001'
    self.archive        = 'elpa-'+self.version+'.tar.gz'
    self.url            = 'https://elpa.mpcdf.mpg.de/software/tarball-archive/Releases/'+self.version+'/'+self.archive
    self.supportssingle = True
    self.fortran        = True
    self.ProcessArgs(argdb)


  def Precondition(self,slepc,petsc):
    pkg = self.packagename.upper()
    if not 'scalapack' in petsc.packages:
      self.log.Exit('The ELPA interface requires that PETSc has been built with ScaLAPACK')
    if petsc.language == 'c++':
      self.log.Exit('The ELPA interface currently does not support compilation with C++')
    package.Package.Precondition(self,slepc,petsc)


  def SampleCode(self,petsc):
    code = '#include <stdlib.h>\n'
    code += '#include <elpa/elpa.h>\n'
    code += 'int main() {\n'
    code += '  elpa_t handle;\n'
    code += '  int error;\n'
    code += '  if (elpa_init(20200417)!=ELPA_OK) exit(1);\n'
    code += '  handle = elpa_allocate(&error);\n'
    code += '  elpa_set(handle,"nev",10,&error);\n'
    code += '  elpa_deallocate(handle,&error);\n'
    code += '  elpa_uninit(&error);\n'
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
      dirs = self.GenerateGuesses('elpa',archdir) + self.GenerateGuesses('elpa',archdir,'lib64')
      incdirs = self.GenerateGuesses('elpa',archdir,'include') + self.GenerateGuesses('elpa',archdir,os.path.join('include','elpa-'+self.version))

    libs = [self.packagelibs] if self.packagelibs else ['-lelpa']
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
          slepcconf.write('#define SLEPC_HAVE_ELPA 1\n')
          slepcvars.write('ELPA_LIB = ' + self.libflags + '\n')
          slepcvars.write('ELPA_INCLUDE = ' + self.includeflags + '\n')
          self.havepackage = True
          self.packageflags = ' '.join(l+f)
          return

    self.log.Exit('Unable to link with ELPA library in directories'+' '.join(dirs)+' with libraries and link flags '+' '.join(libs))


  def DownloadAndInstall(self,slepcconf,slepcvars,slepc,petsc,archdir,prefixdir):
    externdir = slepc.GetExternalPackagesDir(archdir)
    builddir  = self.Download(externdir,slepc.downloaddir)
    incdir,libdir = slepc.CreatePrefixDirs(prefixdir)

    # Check for autoreconf
    (result,output) = self.RunCommand('autoreconf --help')
    if result:
      self.log.Exit('--download-elpa requires that the command autoreconf is available on your PATH')

    # Build package
    confopt = ['--prefix='+prefixdir, '--libdir='+os.path.join(prefixdir,'lib'), 'CC="'+petsc.cc+'"', 'CFLAGS="'+petsc.getCFlags()+'"', 'F77="'+petsc.fc+'"', 'FFLAGS="'+petsc.getFFlags()+'"', 'FC="'+petsc.fc+'"', 'FCFLAGS="'+petsc.getFFlags()+'"', 'CXX="'+petsc.cxx+'"', 'CXXFLAGS="'+petsc.getCXXFlags()+'"', 'CPP="'+petsc.cpp+'"', 'SCALAPACK_LDFLAGS="'+petsc.scalapack_lib+'"', '--disable-sse', '--disable-sse-assembly', '--disable-avx', '--disable-avx2', '--disable-avx512', '--disable-c-tests', '--disable-cpp-tests']
    if petsc.fc_version.startswith('nvf'):
      confopt.append('LIBS="'+petsc.blaslapack_lib+' -lnvf"')
    else:
      confopt.append('LIBS="'+petsc.blaslapack_lib+'"')
    if petsc.mpiuni or petsc.msmpi:
      confopt.append('--with-mpi=no')
    if petsc.precision == 'single':
      confopt.append('--enable-single-precision')
    (result,output) = self.RunCommand('cd '+builddir+'&& ./configure '+' '.join(confopt)+' '+self.buildflags+' && '+petsc.make+' -j'+petsc.make_np+' && '+petsc.make+' install')
    if result:
      self.log.Exit('Installation of ELPA failed')

    # Check build
    code = self.SampleCode(petsc)
    altlibdir = os.path.join(prefixdir,'lib64')
    for ldir in [libdir,altlibdir]:
      if petsc.buildsharedlib:
        l = self.slflag + ldir + ' -L' + ldir + ' -lelpa'
      else:
        l = '-L' + ldir + ' -lelpa'
      f = '-I' + os.path.join(incdir,self.GetDirectoryName())
      (result, output) = self.Link([],[],l+' '+f,code,f,petsc.language)
      if result: break

    if not result:
      self.log.Exit('Unable to link with downloaded ELPA')

    # Write configuration files
    self.libflags = l
    self.includeflags = f
    slepcconf.write('#define SLEPC_HAVE_ELPA 1\n')
    slepcvars.write('ELPA_LIB = ' + self.libflags + '\n')
    slepcvars.write('ELPA_INCLUDE = ' + self.includeflags + '\n')

    self.havepackage = True
    self.packageflags = l+' '+f

