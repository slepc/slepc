/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <petsc/private/ftnimpl.h>
#include <slepcpep.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define pepstoargetinertias_        PEPSTOARGETINERTIAS
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define pepstoargetinertias_        pepstoargetinertias
#endif

SLEPC_EXTERN void pepstoargetinertias_(PEP *pep,PetscInt *nshift,PetscReal *shifts,PetscInt *inertias,PetscErrorCode *ierr)
{
  PetscReal *oshifts;
  PetscInt  *oinertias;
  PetscInt  n;

  CHKFORTRANNULLREAL(shifts);
  CHKFORTRANNULLINTEGER(inertias);
  *ierr = PEPSTOARGetInertias(*pep,&n,&oshifts,&oinertias); if (*ierr) return;
  if (shifts) { *ierr = PetscArraycpy(shifts,oshifts,n); if (*ierr) return; }
  if (inertias) { *ierr = PetscArraycpy(inertias,oinertias,n); if (*ierr) return; }
  *nshift = n;
  *ierr = PetscFree(oshifts);if (*ierr) return;
  *ierr = PetscFree(oinertias);
}
