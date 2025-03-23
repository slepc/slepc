!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Used by slepcsysmod.F90 to create Fortran module file
!

#include "slepc/finclude/slepcsys.h"

      type tSlepcSC
        PetscFortranAddr:: v PETSC_FORTRAN_TYPE_INITIALIZE
      end type tSlepcSC

! Default tolerance for the different solvers, depending on the precision

#if defined(PETSC_USE_REAL_SINGLE)
      PetscReal, parameter :: SLEPC_DEFAULT_TOL =  1e-5
#elif defined(PETSC_USE_REAL_DOUBLE)
      PetscReal, parameter :: SLEPC_DEFAULT_TOL =  1e-8
#elif defined(PETSC_USE_REAL___FLOAT128)
      PetscReal, parameter :: SLEPC_DEFAULT_TOL = 1e-16
#elif defined(PETSC_USE_REAL___FP16)
      PetscReal, parameter :: SLEPC_DEFAULT_TOL =  1e-2
#endif
