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

class Ksvd(package.Package):

  def __init__(self,argdb,log):
    package.Package.__init__(self,argdb,log)
    self.packagename    = 'ksvd'
    self.packagetype    = 'cmake'
    self.installable    = True
    self.downloadable   = True
    self.version        = '2.0.0'
    self.url            = 'https://github.com/ecrc/ksvd/archive/v'+self.version+'.tar.gz'
    self.archive        = 'ksvd-'+self.version+'.tar.gz'
    self.supportsscalar = ['real']
    self.ProcessArgs(argdb)


  def Precondition(self,slepc,petsc):
    self.elpa  = self.Require('elpa')
    self.polar = self.Require('polar')
    if getattr(self,'download',False) and not hasattr(petsc,'cmake'):
      self.log.Exit('The KSVD interface requires CMake for building')
    package.Package.Precondition(self,slepc,petsc)


  def SampleCode(self,petsc):
    code =  '#include <stdlib.h>\n'
    code += '#include "ksvd.h"\n'
    code += 'int main() {\n'
    code += '  int n,i1=1,info,lwork=-1,*w2,liwork=-1,descA[9],descU[9],descVT[9];\n'
    code += '  double *A,*S,*U,*VT,*w1;\n'
    code += '  pdgeqsvd("V","V","r",n,n,A,i1,i1,descA,S,U,i1,i1,descU,VT,i1,i1,descVT,w1,lwork,w2,liwork,&info);\n'
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
      dirs = self.GenerateGuesses('ksvd',archdir) + self.GenerateGuesses('ksvd',archdir,'lib64')
      incdirs = self.GenerateGuesses('ksvd',archdir,'include')

    libs = [self.packagelibs] if self.packagelibs else ['-lksvd']
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
        (result, output) = self.Link([],[],' '.join(l+f+self.elpa.libflags+self.polar.libflags),code,' '.join(f+[self.elpa.includeflags]),petsc.language)
        if result:
          self.libflags = ' '.join(l)
          self.includeflags = ' '.join(f)
          slepcconf.write('#define SLEPC_HAVE_KSVD 1\n')
          slepcvars.write('KSVD_LIB = ' + self.libflags + '\n')
          slepcvars.write('KSVD_INCLUDE = ' + self.includeflags + '\n')
          self.havepackage = True
          self.packageflags = ' '.join(l+f)
          return

    self.log.Exit('Unable to link with KSVD library in directories'+' '.join(dirs)+' with libraries and link flags '+' '.join(libs))


  def DownloadAndInstall(self,slepcconf,slepcvars,slepc,petsc,archdir,prefixdir):
    externdir = slepc.GetExternalPackagesDir(archdir)
    builddir  = self.Download(externdir,slepc.downloaddir)
    incdir,libdir = slepc.CreatePrefixDirs(prefixdir)

    # Download CMake module files
    (result, output) = self.RunCommand('cd '+os.path.join(builddir,'cmake_modules')+' && rm -rf ecrc && git clone https://github.com/ecrc/ecrc_cmake.git ecrc')
    if result:
      self.log.Exit('Unable to download CMake module files needed for KSVD')

    # Patch FindELPA.cmake CMake module file
    (result,output) = self.RunCommand('cd '+os.path.join(builddir,'cmake_modules','ecrc','modules','find')+' && '+petsc.sedinplace+' '+'-e "s?NAMES elpa.h?NAMES elpa/elpa.h?" -e "s?elpa_dgetrf?elpa_init?" FindELPA.cmake')
    if result:
      self.log.Exit('Problem when patching file FindELPA.cmake')

    # Patch pdgeqsvd.c to use the API of recent ELPA
    fname = os.path.join(builddir,'src','pdgeqsvd.c')
    oldcode1 = '''int useQr, THIS_REAL_ELPA_KERNEL_API;
        int mpi_comm_rows, mpi_comm_cols;
        int mpierr = elpa_get_communicators(MPI_Comm_c2f(MPI_COMM_WORLD), myrow, mycol, &mpi_comm_rows, &mpi_comm_cols);
        useQr = 0;
        THIS_REAL_ELPA_KERNEL_API = ELPA2_REAL_KERNEL_AVX_BLOCK6;
        *info = elpa_solve_evp_real_2stage( n, n, U, mloc, 
                                            S, VT, 
                                            mloc, nb, nloc, 
                                            mpi_comm_rows, mpi_comm_cols, MPI_Comm_c2f(MPI_COMM_WORLD),
                                            THIS_REAL_ELPA_KERNEL_API, useQr);'''
    newcode1 = '''elpa_t handle;
        int error_elpa;
        if (elpa_init(20200417) != ELPA_OK) {
          fprintf(stderr, "Error: ELPA API version not supported");
          exit(1);
        }
        handle = elpa_allocate(&error_elpa);
        elpa_set(handle,"na",n,&error_elpa); assert(error_elpa == ELPA_OK);
        elpa_set(handle,"nev",n,&error_elpa); assert(error_elpa == ELPA_OK);
        elpa_set(handle,"local_nrows",mloc,&error_elpa); assert(error_elpa == ELPA_OK);
        elpa_set(handle,"local_ncols",nloc,&error_elpa); assert(error_elpa == ELPA_OK);
        elpa_set(handle,"nblk",nb,&error_elpa); assert(error_elpa == ELPA_OK);
        elpa_set(handle,"mpi_comm_parent",MPI_Comm_c2f(MPI_COMM_WORLD),&error_elpa); assert(error_elpa == ELPA_OK);
        elpa_set(handle,"process_row",myrow,&error_elpa); assert(error_elpa == ELPA_OK);
        elpa_set(handle,"process_col",mycol,&error_elpa); assert(error_elpa == ELPA_OK);
        error_elpa = elpa_setup(handle); assert(error_elpa == ELPA_OK);
        elpa_set(handle,"solver",ELPA_SOLVER_2STAGE,&error_elpa); assert(error_elpa == ELPA_OK);
        elpa_eigenvectors(handle,U,S,VT,&error_elpa);  assert(error_elpa == ELPA_OK);
        elpa_deallocate(handle, &error_elpa); assert(error_elpa == ELPA_OK);
        elpa_uninit(&error_elpa); assert(error_elpa == ELPA_OK);'''
    oldcode2 = r'''#include "ksvd.h"'''
    newcode2 = r'''#include "ksvd.h"
int pdgeqdwh(char *jobh,int m,int n,double *A,int iA,int jA,int *descA,double *H,int iH,int jH,int *descH,double *Work1,int lWork1,double *Work2,int lWork2,int *info);
int pdgezolopd(char *jobh,int m,int n,double *A,int iA,int jA,int *descA,double *H,int iH,int jH,int *descH,double *Work1,int lWork1,double *Work2,int lWork2,int *info);'''
    with open(fname,'r') as file:
      sourcecode = file.read()
    sourcecode = sourcecode.replace(oldcode1,newcode1).replace(oldcode2,newcode2)
    with open(fname,'w') as file:
      file.write(sourcecode)

    # Build package
    builddir = slepc.CreateDir(builddir,'build')
    confopt = ['-DCMAKE_INSTALL_PREFIX='+prefixdir, '-DCMAKE_INSTALL_NAME_DIR:STRING="'+os.path.join(prefixdir,'lib')+'"', '-DCMAKE_C_COMPILER="'+petsc.cc+'"', '-DCMAKE_C_FLAGS:STRING="'+petsc.getCFlags()+'"', '-DELPA_INCDIR="'+os.path.join(incdir,'elpa-'+self.elpa.version)+'"', '-DELPA_LIBDIR="'+libdir+'"', '-DPOLAR_DIR="'+prefixdir+'"', '-DBLAS_LIBRARIES="'+petsc.blaslapack_lib+'"']
    confopt.append('-DCMAKE_BUILD_TYPE='+('Debug' if petsc.debug else 'Release'))
    if petsc.buildsharedlib:
      confopt = confopt + ['-DBUILD_SHARED_LIBS=ON', '-DCMAKE_INSTALL_RPATH:PATH='+os.path.join(prefixdir,'lib')]
    else:
      confopt.append('-DBUILD_SHARED_LIBS=OFF')
    if 'MSYSTEM' in os.environ:
      confopt.append('-G "MSYS Makefiles"')
    (result,output) = self.RunCommand('cd '+builddir+' && '+petsc.cmake+' '+' '.join(confopt)+' '+self.buildflags+' .. && '+petsc.make+' -j'+petsc.make_np+' && '+petsc.make+' install')
    if result:
      self.log.Exit('Installation of KSVD failed')

    # Patch include file
    (result,output) = self.RunCommand('cd '+incdir+' && '+petsc.sedinplace+' '+'-e "/myscalapack.h/d" -e "/flops.h/d" ksvd.h')
    if result:
      self.log.Exit('Problem when patching include file ksvd.h')
    fname = os.path.join(incdir,'ksvd.h')
    oldcode1 = r'''int pdgeqsvd( char *jobu, char *jobvt, char *eigtype, '''
    newcode1 = r'''int pdgeqsvd( const char *jobu, const char *jobvt, const char *eigtype,'''
    with open(fname,'r') as file:
      sourcecode = file.read()
    sourcecode = sourcecode.replace(oldcode1,newcode1)
    with open(fname,'w') as file:
      file.write(sourcecode)

    # Check build
    code = self.SampleCode(petsc)
    altlibdir = os.path.join(prefixdir,'lib64')
    for ldir in [libdir,altlibdir]:
      if petsc.buildsharedlib:
        l = self.slflag + ldir + ' -L' + ldir + ' -lksvd'
      else:
        l = '-L' + ldir + ' -lksvd'
      f = '-I' + incdir
      (result, output) = self.Link([],[],l+' '+f+' '+self.elpa.libflags+' '+self.polar.libflags,code,f+' '+self.elpa.includeflags,petsc.language)
      if result: break

    if not result:
      self.log.Exit('Unable to link with downloaded KSVD')

    # Write configuration files
    self.libflags = l
    self.includeflags = f
    slepcconf.write('#define SLEPC_HAVE_KSVD 1\n')
    slepcvars.write('KSVD_LIB = ' + self.libflags + '\n')
    slepcvars.write('KSVD_INCLUDE = ' + self.includeflags + '\n')

    self.havepackage = True
    self.packageflags = l+' '+f

