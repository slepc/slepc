#
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#  SLEPc - Scalable Library for Eigenvalue Problem Computations
#  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
#
#  This file is part of SLEPc.
#  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
#  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#

import log, package

class Feast(package.Package):

  def __init__(self,argdb,log,petscpackages):
    package.Package.__init__(self,argdb,log)
    self.packagename    = 'feast'
    self.installable    = True
    self.petscdepend    = 'mkl'
    self.supports64bint = True
    self.supportsprecis.extend(['single','__float128'])
    self.ProcessArgs(argdb,petscpackages)

  def SampleCode(self,petsc):
    if petsc.scalar == 'real':
      if petsc.precision == 'single':
        function = 'sfeast_srci'
        rtype = 'float'
        ctype = 'MKL_Complex8'
      else:
        function = 'dfeast_srci'
        rtype = 'double'
        ctype = 'MKL_Complex16'
      stype = rtype
    else:
      if petsc.precision == 'single':
        function = 'cfeast_hrci'
        rtype = 'float'
        ctype = 'MKL_Complex8'
      else:
        function = 'zfeast_hrci'
        rtype = 'double'
        ctype = 'MKL_Complex16'
      stype = ctype

    code = '#include <mkl.h>\n'
    code += 'int main() {\n'
    code += '  ' + rtype + ' epsout=0.0,*evals=NULL,*errest=NULL,inta,intb;\n'
    code += '  ' + ctype + ' Ze,*work2=NULL;\n'
    code += '  ' + stype + ' *Aq=NULL,*Bq=NULL,*pV=NULL,*work1=NULL;\n'
    code += '  MKL_INT fpm[128],ijob,n,loop,ncv,nconv,info;\n'
    code += '  feastinit(fpm);\n'
    code += '  ' + function + '(&ijob,&n,&Ze,work1,work2,Aq,Bq,fpm,&epsout,&loop,&inta,&intb,&ncv,evals,pV,&nconv,errest,&info);\n'
    code += '  return 0;\n}\n'
    return code

  def Check(self,slepcconf,slepcvars,petsc,archdir):
    if not 'mkl' in petsc.packages:
      self.log.Exit('The FEAST interface requires that PETSc has been built with Intel MKL (libraries and includes)')
    if petsc.precision == '__float128':
      self.skippackage = True
      self.log.write('Disabling FEAST interface due to precision = __float128')
    else:
      code = self.SampleCode(petsc)

      (result,output) = self.Link([],[],'',code)
      if not result:
        self.log.write('WARNING: Unable to link with FEAST, maybe your MKL version does not contain it')
        self.log.write('If you do not want to check for FEAST, rerun configure adding --with-feast=0')
        self.havepackage = False
      else:
        slepcconf.write('#define SLEPC_HAVE_FEAST 1\n')
        self.havepackage = True
