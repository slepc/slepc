!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Used by slepcsvdmod.F90 to create Fortran module file
!
#include "slepc/finclude/slepcsvd.h"

      type, extends(tPetscObject) :: tSVD
      end type tSVD

      SVD, parameter :: SLEPC_NULL_SVD = tSVD(0)

      PetscEnum, parameter :: SVD_CONVERGED_TOL          =  1
      PetscEnum, parameter :: SVD_CONVERGED_USER         =  2
      PetscEnum, parameter :: SVD_CONVERGED_MAXIT        =  3
      PetscEnum, parameter :: SVD_DIVERGED_ITS           = -1
      PetscEnum, parameter :: SVD_DIVERGED_BREAKDOWN     = -2
      PetscEnum, parameter :: SVD_CONVERGED_ITERATING    =  0

      PetscEnum, parameter :: SVD_STANDARD               =  1
      PetscEnum, parameter :: SVD_GENERALIZED            =  2
      PetscEnum, parameter :: SVD_HYPERBOLIC             =  3

      PetscEnum, parameter :: SVD_LARGEST                =  0
      PetscEnum, parameter :: SVD_SMALLEST               =  1

      PetscEnum, parameter :: SVD_ERROR_ABSOLUTE         =  0
      PetscEnum, parameter :: SVD_ERROR_RELATIVE         =  1
      PetscEnum, parameter :: SVD_ERROR_NORM             =  2

      PetscEnum, parameter :: SVD_CONV_ABS               =  0
      PetscEnum, parameter :: SVD_CONV_REL               =  1
      PetscEnum, parameter :: SVD_CONV_NORM              =  2
      PetscEnum, parameter :: SVD_CONV_MAXIT             =  3
      PetscEnum, parameter :: SVD_CONV_USER              =  4

      PetscEnum, parameter :: SVD_STOP_BASIC             =  0
      PetscEnum, parameter :: SVD_STOP_USER              =  1
      PetscEnum, parameter :: SVD_STOP_THRESHOLD         =  2

      PetscEnum, parameter :: SVD_TRLANCZOS_GBIDIAG_SINGLE =  0
      PetscEnum, parameter :: SVD_TRLANCZOS_GBIDIAG_UPPER  =  1
      PetscEnum, parameter :: SVD_TRLANCZOS_GBIDIAG_LOWER  =  2

      PetscEnum, parameter :: SVD_PRIMME_HYBRID          =  1
      PetscEnum, parameter :: SVD_PRIMME_NORMALEQUATIONS =  2
      PetscEnum, parameter :: SVD_PRIMME_AUGMENTED       =  3

!
!   Possible arguments to SVDMonitorSet()
!
      external SVDMONITORFIRST
      external SVDMONITORALL
      external SVDMONITORCONVERGED

      external SVDMonitorConvergedDestroy

#if defined(_WIN32) && defined(PETSC_USE_SHARED_LIBRARIES)
!DEC$ ATTRIBUTES DLLEXPORT::SLEPC_NULL_SVD
#endif
