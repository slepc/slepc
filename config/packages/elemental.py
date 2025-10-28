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

class Elemental(package.Package):

  def __init__(self,argdb,log,petscpackages):
    package.Package.__init__(self,argdb,log)
    self.packagename    = 'elemental'
    self.installable    = True
    self.petscdepend    = 'elemental'
    self.supports64bint = True
    self.supportsprecis.append('single')
    self.ProcessArgs(argdb,petscpackages)

  def Check(self,slepcconf,slepcvars,petsc,archdir):
    if not 'elemental' in petsc.packages:
      self.log.Exit('The Elemental interface requires that PETSc has been built with Elemental')

    code =  '#include <petscmatelemental.h>\n'
    code += 'int main() {\n'
    code += '  El::mpi::Comm comm = El::mpi::COMM_WORLD;\n'
    code += '  El::Grid grid( comm );\n'
    code += '  El::DistMatrix<PetscReal,El::VR,El::STAR> w( grid );\n'
    code += '  El::DistMatrix<PetscElemScalar> H( 10, 10, grid );\n'
    code += '  El::HermitianEig( El::LOWER, H, w, H );\n'
    code += '  return 0;\n}\n'

    (result,output) = self.Link([],[],'',code,clanguage='c++')
    if not result:
      self.log.write('WARNING: Unable to link with Elemental')
      self.log.write('If you do not want to check for Elemental, rerun configure adding --with-elemental=0')
      self.havepackage = False
    else:
      slepcconf.write('#define SLEPC_HAVE_ELEMENTAL 1\n')
      self.havepackage = True
