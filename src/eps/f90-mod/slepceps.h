!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!  SLEPc is distributed under a 2-clause BSD license (see LICENSE).
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!  Used by slepcepsmod.F90 to create Fortran module file
!
#include "slepc/finclude/slepceps.h"

      type, extends(tPetscObject) :: tEPS
      end type tEPS

      EPS, parameter :: SLEPC_NULL_EPS = tEPS(0)

      PetscEnum, parameter :: EPS_CONVERGED_TOL          =  1
      PetscEnum, parameter :: EPS_CONVERGED_USER         =  2
      PetscEnum, parameter :: EPS_DIVERGED_ITS           = -1
      PetscEnum, parameter :: EPS_DIVERGED_BREAKDOWN     = -2
      PetscEnum, parameter :: EPS_DIVERGED_SYMMETRY_LOST = -3
      PetscEnum, parameter :: EPS_CONVERGED_ITERATING    =  0

      PetscEnum, parameter :: EPS_HEP                    =  1
      PetscEnum, parameter :: EPS_GHEP                   =  2
      PetscEnum, parameter :: EPS_NHEP                   =  3
      PetscEnum, parameter :: EPS_GNHEP                  =  4
      PetscEnum, parameter :: EPS_PGNHEP                 =  5
      PetscEnum, parameter :: EPS_GHIEP                  =  6
      PetscEnum, parameter :: EPS_BSE                    =  7

      PetscEnum, parameter :: EPS_LARGEST_MAGNITUDE      =  1
      PetscEnum, parameter :: EPS_SMALLEST_MAGNITUDE     =  2
      PetscEnum, parameter :: EPS_LARGEST_REAL           =  3
      PetscEnum, parameter :: EPS_SMALLEST_REAL          =  4
      PetscEnum, parameter :: EPS_LARGEST_IMAGINARY      =  5
      PetscEnum, parameter :: EPS_SMALLEST_IMAGINARY     =  6
      PetscEnum, parameter :: EPS_TARGET_MAGNITUDE       =  7
      PetscEnum, parameter :: EPS_TARGET_REAL            =  8
      PetscEnum, parameter :: EPS_TARGET_IMAGINARY       =  9
      PetscEnum, parameter :: EPS_ALL                    = 10
      PetscEnum, parameter :: EPS_WHICH_USER             = 11

      PetscEnum, parameter :: EPS_BALANCE_NONE           =  0
      PetscEnum, parameter :: EPS_BALANCE_ONESIDE        =  1
      PetscEnum, parameter :: EPS_BALANCE_TWOSIDE        =  2
      PetscEnum, parameter :: EPS_BALANCE_USER           =  3

      PetscEnum, parameter :: EPS_RITZ                   =  0
      PetscEnum, parameter :: EPS_HARMONIC               =  1
      PetscEnum, parameter :: EPS_HARMONIC_RELATIVE      =  2
      PetscEnum, parameter :: EPS_HARMONIC_RIGHT         =  3
      PetscEnum, parameter :: EPS_HARMONIC_LARGEST       =  4
      PetscEnum, parameter :: EPS_REFINED                =  5
      PetscEnum, parameter :: EPS_REFINED_HARMONIC       =  6

      PetscEnum, parameter :: EPS_ERROR_ABSOLUTE         =  0
      PetscEnum, parameter :: EPS_ERROR_RELATIVE         =  1
      PetscEnum, parameter :: EPS_ERROR_BACKWARD         =  2

      PetscEnum, parameter :: EPS_CONV_ABS               =  0
      PetscEnum, parameter :: EPS_CONV_REL               =  1
      PetscEnum, parameter :: EPS_CONV_NORM              =  2
      PetscEnum, parameter :: EPS_CONV_USER              =  3

      PetscEnum, parameter :: EPS_STOP_BASIC             =  0
      PetscEnum, parameter :: EPS_STOP_USER              =  1
      PetscEnum, parameter :: EPS_STOP_THRESHOLD         =  2

      PetscEnum, parameter :: EPS_POWER_SHIFT_CONSTANT   =  0
      PetscEnum, parameter :: EPS_POWER_SHIFT_RAYLEIGH   =  1
      PetscEnum, parameter :: EPS_POWER_SHIFT_WILKINSON  =  2

      PetscEnum, parameter :: EPS_KRYLOVSCHUR_BSE_SHAO         =  0
      PetscEnum, parameter :: EPS_KRYLOVSCHUR_BSE_GRUNING      =  1
      PetscEnum, parameter :: EPS_KRYLOVSCHUR_BSE_PROJECTEDBSE =  2

      PetscEnum, parameter :: EPS_LANCZOS_REORTHOG_LOCAL     =  0
      PetscEnum, parameter :: EPS_LANCZOS_REORTHOG_FULL      =  1
      PetscEnum, parameter :: EPS_LANCZOS_REORTHOG_SELECTIVE =  2
      PetscEnum, parameter :: EPS_LANCZOS_REORTHOG_PERIODIC  =  3
      PetscEnum, parameter :: EPS_LANCZOS_REORTHOG_PARTIAL   =  4
      PetscEnum, parameter :: EPS_LANCZOS_REORTHOG_DELAYED   =  5

      PetscEnum, parameter :: EPS_PRIMME_DYNAMIC             =  1
      PetscEnum, parameter :: EPS_PRIMME_DEFAULT_MIN_TIME    =  2
      PetscEnum, parameter :: EPS_PRIMME_DEFAULT_MIN_MATVECS =  3
      PetscEnum, parameter :: EPS_PRIMME_ARNOLDI             =  4
      PetscEnum, parameter :: EPS_PRIMME_GD                  =  5
      PetscEnum, parameter :: EPS_PRIMME_GD_PLUSK            =  6
      PetscEnum, parameter :: EPS_PRIMME_GD_OLSEN_PLUSK      =  7
      PetscEnum, parameter :: EPS_PRIMME_JD_OLSEN_PLUSK      =  8
      PetscEnum, parameter :: EPS_PRIMME_RQI                 =  9
      PetscEnum, parameter :: EPS_PRIMME_JDQR                = 10
      PetscEnum, parameter :: EPS_PRIMME_JDQMR               = 11
      PetscEnum, parameter :: EPS_PRIMME_JDQMR_ETOL          = 12
      PetscEnum, parameter :: EPS_PRIMME_SUBSPACE_ITERATION  = 13
      PetscEnum, parameter :: EPS_PRIMME_LOBPCG_ORTHOBASIS   = 14
      PetscEnum, parameter :: EPS_PRIMME_LOBPCG_ORTHOBASISW  = 15

      PetscEnum, parameter :: EPS_CISS_QUADRULE_TRAPEZOIDAL  =  1
      PetscEnum, parameter :: EPS_CISS_QUADRULE_CHEBYSHEV    =  2

      PetscEnum, parameter :: EPS_CISS_EXTRACTION_RITZ       =  0
      PetscEnum, parameter :: EPS_CISS_EXTRACTION_HANKEL     =  1

      PetscEnum, parameter :: EPS_EVSL_DOS_KPM               =  0
      PetscEnum, parameter :: EPS_EVSL_DOS_LANCZOS           =  1

      PetscEnum, parameter :: EPS_EVSL_DAMPING_NONE          =  0
      PetscEnum, parameter :: EPS_EVSL_DAMPING_JACKSON       =  1
      PetscEnum, parameter :: EPS_EVSL_DAMPING_SIGMA         =  2

!
!   Possible arguments to EPSMonitorSet()
!
      external EPSMONITORFIRST
      external EPSMONITORALL
      external EPSMONITORCONVERGED

      external EPSMonitorConvergedDestroy

#if defined(_WIN32) && defined(PETSC_USE_SHARED_LIBRARIES)
!DEC$ ATTRIBUTES DLLEXPORT::SLEPC_NULL_EPS
#endif
