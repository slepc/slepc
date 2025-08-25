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

class Polar(package.Package):

  def __init__(self,argdb,log):
    package.Package.__init__(self,argdb,log)
    self.packagename    = 'polar'
    self.packagetype    = 'cmake'
    self.installable    = True
    self.downloadable   = True
    self.commit         = 'ae12cb521b2174ba31fa9d5adadbacb85819381a'
    self.url            = 'https://github.com/ecrc/polar/archive/'+self.commit+'.tar.gz'
    self.archive        = 'polar-'+self.commit+'.tar.gz'
    self.supportsscalar = ['real']
    self.ProcessArgs(argdb)


  def Precondition(self,slepc,petsc):
    pkg = self.packagename.upper()
    if not 'mkl' in petsc.packages:
      self.log.Exit('The '+pkg+' interface requires that PETSc has been built with Intel MKL (libraries and includes)')
    if getattr(self,'download',False) and not hasattr(petsc,'cmake'):
      self.log.Exit('The POLAR interface requires CMake for building')
    package.Package.Precondition(self,slepc,petsc)


  def SampleCode(self,petsc):
    code =  '#include <stdlib.h>\n'
    code += '#include "polar.h"\n'
    code += 'int main() {\n'
    code += '  int n,i1=1,info,lwork=-1,descA[9],descH[9];\n'
    code += '  double *A,*H,*w1,*w2;\n'
    code += '  pdgeqdwh("H",n,n,A,i1,i1,descA,H,i1,i1,descH,w1,lwork,w2,lwork,&info);\n'
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
      dirs = self.GenerateGuesses('polar',archdir) + self.GenerateGuesses('polar',archdir,'lib64')
      incdirs = self.GenerateGuesses('polar',archdir,'include')

    libs = [self.packagelibs] if self.packagelibs else ['-lpolar']
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
          slepcconf.write('#define SLEPC_HAVE_POLAR 1\n')
          slepcvars.write('POLAR_LIB = ' + self.libflags + '\n')
          slepcvars.write('POLAR_INCLUDE = ' + self.includeflags + '\n')
          self.havepackage = True
          self.packageflags = ' '.join(l+f)
          return

    self.log.Exit('Unable to link with POLAR library in directories'+' '.join(dirs)+' with libraries and link flags '+' '.join(libs))


  def DownloadAndInstall(self,slepcconf,slepcvars,slepc,petsc,archdir,prefixdir):
    externdir = slepc.GetExternalPackagesDir(archdir)
    builddir  = self.Download(externdir,slepc.downloaddir)
    incdir,libdir = slepc.CreatePrefixDirs(prefixdir)

    # Download CMake module files
    (result, output) = self.RunCommand('cd '+os.path.join(builddir,'cmake_modules')+' && rm -rf ecrc && git clone https://github.com/ecrc/ecrc_cmake.git ecrc')
    if result:
      self.log.Exit('Unable to download CMake module files needed for POLAR')

    # Patch pdgeqdwh.c to avoid output even with verbose=false
    fname = os.path.join(builddir,'src','pdgeqdwh.c')
    oldcode1 = r'''if (myrank_mpi == 0) { fprintf(stderr, "\nItConv %d itcqr %d itcpo %d norm_est %2.4e Li %2.4e \n", itconv, itcqr, itcpo, norm_est, Li); fprintf(stderr, "It Facto Conv\n");}'''
    newcode1 = r'''if (verbose && myrank_mpi == 0) { fprintf(stderr, "\nItConv %d itcqr %d itcpo %d norm_est %2.4e Li %2.4e \n", itconv, itcqr, itcpo, norm_est, Li); fprintf(stderr, "It Facto Conv\n");}'''
    oldcode2 = r'''if (myrank_mpi == 0) {
        fprintf(stderr, "#\n");'''
    newcode2 = r'''if (verbose && myrank_mpi == 0) {
        fprintf(stderr, "#\n");'''
    with open(fname,'r') as file:
      sourcecode = file.read()
    sourcecode = sourcecode.replace(oldcode1,newcode1).replace(oldcode2,newcode2)
    with open(fname,'w') as file:
      file.write(sourcecode)

    # Patch pdgezolopd.c to avoid output even with verbose=false
    fname = os.path.join(builddir,'src','pdgezolopd.c')
    oldcode1 = '''if ( myrank_mpi == 0 ){ 
        fprintf(stderr, " The number of subproblems to be solved independently is %d'''
    newcode1 = '''if (verbose && myrank_mpi == 0 ){ 
        fprintf(stderr, " The number of subproblems to be solved independently is %d'''
    oldcode2 = r'''if ( myrank_mpi == 0 ) {
        fprintf(stderr, "#\n");'''
    newcode2 = r'''if (verbose && myrank_mpi == 0 ) {
        fprintf(stderr, "#\n");'''
    oldcode3 = r'''        pxerbla_( ictxt_all, "PDGEZOLOPD", &(int){-1*info[0]} ); '''
    newcode3 = r'''        pxerbla_( &ictxt_all, "PDGEZOLOPD", &(int){-1*info[0]} ); '''
    with open(fname,'r') as file:
      sourcecode = file.read()
    sourcecode = sourcecode.replace(oldcode1,newcode1).replace(oldcode2,newcode2).replace(oldcode3,newcode3)
    with open(fname,'w') as file:
      file.write(sourcecode)

    # Patch include files
    fname = os.path.join(builddir,'include','myscalapack.h')
    oldcode1 = '''#define descinit_ descinit'''
    newcode1 = '''#define descinit_ descinit
#define pdgemr2d_ pdgemr2d
'''
    oldcode2 = r'''extern int numroc_( int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);'''
    newcode2 = r'''extern int numroc_( int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
extern void pdgemr2d_(int *m, int *n, double *a, int *ia, int *ja, int *desca, double *b, int *ib, int *jb, int *descb, int *ictxt );
extern void Cblacs_gridmap(int *ConTxt, int *usermap, int ldup, int nprow0, int npcol0);
extern int Csys2blacs_handle(MPI_Comm SysCtxt);
'''
    with open(fname,'r') as file:
      sourcecode = file.read()
    sourcecode = sourcecode.replace(oldcode1,newcode1).replace(oldcode2,newcode2)
    with open(fname,'w') as file:
      file.write(sourcecode)

    fname = os.path.join(builddir,'include','polar.h')
    oldcode1 = r'''int mellipke( double alpha,  
              double *k, double *e);'''
    newcode1 = r'''int mellipke( double alpha,  
              double *k, double *e);
int choosem(double con, int *m);'''
    with open(fname,'r') as file:
      sourcecode = file.read()
    sourcecode = sourcecode.replace(oldcode1,newcode1)
    with open(fname,'w') as file:
      file.write(sourcecode)

    # Build package
    builddir = slepc.CreateDir(builddir,'build')
    confopt = ['-DCMAKE_INSTALL_PREFIX='+prefixdir, '-DCMAKE_INSTALL_NAME_DIR:STRING="'+os.path.join(prefixdir,'lib')+'"', '-DCMAKE_C_COMPILER="'+petsc.cc+'"', '-DCMAKE_C_FLAGS:STRING="'+petsc.getCFlags()+'"']
    confopt.append('-DCMAKE_BUILD_TYPE='+ ('Debug' if petsc.debug else 'Release'))
    if petsc.buildsharedlib:
      confopt = confopt + ['-DBUILD_SHARED_LIBS=ON', '-DCMAKE_INSTALL_RPATH:PATH='+os.path.join(prefixdir,'lib')]
    else:
      confopt.append('-DBUILD_SHARED_LIBS=OFF')
    if 'MSYSTEM' in os.environ:
      confopt.append('-G "MSYS Makefiles"')
    (result,output) = self.RunCommand('cd '+builddir+' && '+petsc.cmake+' '+' '.join(confopt)+' '+self.buildflags+' .. && '+petsc.make+' -j'+petsc.make_np+' && '+petsc.make+' install')
    if result:
      self.log.Exit('Installation of POLAR failed')

    # Patch include file
    (result,output) = self.RunCommand('cd '+incdir+' && '+petsc.sedinplace+' '+'-e "/myscalapack.h/d" -e "/flops.h/d" polar.h')
    if result:
      self.log.Exit('Problem when patching include file polar.h')

    # Check build
    code = self.SampleCode(petsc)
    altlibdir = os.path.join(prefixdir,'lib64')
    for ldir in [libdir,altlibdir]:
      if petsc.buildsharedlib:
        l = self.slflag + ldir + ' -L' + ldir + ' -lpolar'
      else:
        l = '-L' + ldir + ' -lpolar'
      f = '-I' + incdir
      (result, output) = self.Link([],[],l+' '+f,code,f,petsc.language)
      if result: break

    if not result:
      self.log.Exit('Unable to link with downloaded POLAR')

    # Write configuration files
    self.libflags = l
    self.includeflags = f
    slepcconf.write('#define SLEPC_HAVE_POLAR 1\n')
    slepcvars.write('POLAR_LIB = ' + self.libflags + '\n')
    slepcvars.write('POLAR_INCLUDE = ' + self.includeflags + '\n')

    self.havepackage = True
    self.packageflags = l+' '+f

