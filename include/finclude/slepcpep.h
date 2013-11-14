!
!  Include file for Fortran use of the PEP object in SLEPc
!
!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SLEPc - Scalable Library for Eigenvalue Problem Computations
!  Copyright (c) 2002-2013, Universitat Politecnica de Valencia, Spain
!
!  This file is part of SLEPc.
!
!  SLEPc is free software: you can redistribute it and/or modify it under  the
!  terms of version 3 of the GNU Lesser General Public License as published by
!  the Free Software Foundation.
!
!  SLEPc  is  distributed in the hope that it will be useful, but WITHOUT  ANY
!  WARRANTY;  without even the implied warranty of MERCHANTABILITY or  FITNESS
!  FOR  A  PARTICULAR PURPOSE. See the GNU Lesser General Public  License  for
!  more details.
!
!  You  should have received a copy of the GNU Lesser General  Public  License
!  along with SLEPc. If not, see <http://www.gnu.org/licenses/>.
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
#include "finclude/slepcpepdef.h"

!  Convergence flags.
!  They should match the flags in $SLEPC_DIR/include/slepcpep.h

      PetscEnum PEP_CONVERGED_TOL
      PetscEnum PEP_DIVERGED_ITS
      PetscEnum PEP_DIVERGED_BREAKDOWN
      PetscEnum PEP_CONVERGED_ITERATING

      parameter (PEP_CONVERGED_TOL          =  2)
      parameter (PEP_DIVERGED_ITS           = -3)
      parameter (PEP_DIVERGED_BREAKDOWN     = -4)
      parameter (PEP_CONVERGED_ITERATING    =  0)

      PetscEnum PEP_GENERAL
      PetscEnum PEP_HERMITIAN
      PetscEnum PEP_GYROSCOPIC

      parameter (PEP_GENERAL                =  1)
      parameter (PEP_HERMITIAN              =  2)
      parameter (PEP_GYROSCOPIC             =  3)

      PetscEnum PEP_LARGEST_MAGNITUDE
      PetscEnum PEP_SMALLEST_MAGNITUDE
      PetscEnum PEP_LARGEST_REAL
      PetscEnum PEP_SMALLEST_REAL
      PetscEnum PEP_LARGEST_IMAGINARY
      PetscEnum PEP_SMALLEST_IMAGINARY
      PetscEnum PEP_TARGET_MAGNITUDE
      PetscEnum PEP_TARGET_REAL
      PetscEnum PEP_TARGET_IMAGINARY

      parameter (PEP_LARGEST_MAGNITUDE      =  1)
      parameter (PEP_SMALLEST_MAGNITUDE     =  2)
      parameter (PEP_LARGEST_REAL           =  3)
      parameter (PEP_SMALLEST_REAL          =  4)
      parameter (PEP_LARGEST_IMAGINARY      =  5)
      parameter (PEP_SMALLEST_IMAGINARY     =  6)
      parameter (PEP_TARGET_MAGNITUDE       =  7)
      parameter (PEP_TARGET_REAL            =  8)
      parameter (PEP_TARGET_IMAGINARY       =  9)

!
!   Possible arguments to PEPMonitorSet()
!
      external PEPMONITORALL
      external PEPMONITORLG
      external PEPMONITORLGALL
      external PEPMONITORCONVERGED
      external PEPMONITORFIRST

!
!  End of Fortran include file for the PEP package in SLEPc
!
