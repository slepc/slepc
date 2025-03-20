/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <petsc/private/ftnimpl.h>
#include <slepcds.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define dspepgetcoefficients_    DSPEPGETCOEFFICIENTS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define dspepgetcoefficients_    dspepgetcoefficients
#endif

SLEPC_EXTERN void dspepgetcoefficients_(DS *ds,PetscReal *pbc,PetscErrorCode *ierr)
{
  PetscReal *opbc;
  PetscInt  d;

  CHKFORTRANNULLREAL(pbc);
  *ierr = DSPEPGetCoefficients(*ds,&opbc); if (*ierr) return;
  *ierr = DSPEPGetDegree(*ds,&d); if (*ierr) return;
  *ierr = PetscArraycpy(pbc,opbc,3*(d+1)); if (*ierr) return;
  *ierr = PetscFree(opbc);
}
