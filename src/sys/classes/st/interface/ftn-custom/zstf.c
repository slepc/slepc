/*
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   SLEPc - Scalable Library for Eigenvalue Problem Computations
   Copyright (c) 2002-, Universitat Politecnica de Valencia, Spain

   This file is part of SLEPc.
   SLEPc is distributed under a 2-clause BSD license (see LICENSE).
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
*/

#include <petsc/private/fortranimpl.h>
#include <slepcst.h>

#if defined(PETSC_HAVE_FORTRAN_CAPS)
#define stdestroy_                STDESTROY
#elif !defined(PETSC_HAVE_FORTRAN_UNDERSCORE)
#define stdestroy_                stdestroy
#endif

SLEPC_EXTERN void stdestroy_(ST *st,PetscErrorCode *ierr)
{
  PETSC_FORTRAN_OBJECT_F_DESTROYED_TO_C_NULL(st);
  *ierr = STDestroy(st); if (*ierr) return;
  PETSC_FORTRAN_OBJECT_C_NULL_TO_F_DESTROYED(st);
}
