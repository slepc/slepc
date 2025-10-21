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

class Scalapack(package.Package):

  def __init__(self,argdb,log,petscpackages):
    package.Package.__init__(self,argdb,log)
    self.packagename    = 'scalapack'
    self.installable    = True
    self.petscdepend    = 'scalapack'
    self.supports64bint = True
    self.supportsprecis.extend(['single','__float128'])
    self.ProcessArgs(argdb,petscpackages)

  def Functions(self,petsc):
    if petsc.scalar == 'real':
      if petsc.precision == 'single':
        functions = ['pssyev','pssygvx','psgesvd']
      else:
        functions = ['pdsyev','pdsygvx','pdgesvd']
    else:
      if petsc.precision == 'single':
        functions = ['pcheev','pchegvx','pcgesvd']
      else:
        functions = ['pzheev','pzhegvx','pzgesvd']
    return functions

  def Check(self,slepcconf,slepcvars,petsc,archdir):
    if not 'scalapack' in petsc.packages:
      self.log.Exit('The ScaLAPACK interface requires that PETSc has been built with ScaLAPACK')

    if petsc.precision == '__float128':
      self.skippackage = True
      self.log.write('Disabling ScaLAPACK interface due to precision = __float128')
    else:
      functions = self.Functions(petsc)
      self.FortranLib(slepcconf,slepcvars,[''],'',functions)

